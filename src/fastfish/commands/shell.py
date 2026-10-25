"""CLASS-inspired command prompt for the shared FAST reduction engine."""
from fastfish._paths import PROJECT_ROOT
from pathlib import Path
import argparse,cmd,fcntl,json,os,shlex,signal,socket,subprocess,sys,time,webbrowser
import numpy as np
from fastfish.processing.session import read_session,write_session
from fastfish.calibration.reference import smoothing_rows
from fastfish.baseline.windows import normalize_windows
from fastfish.baseline.components import default_components,normalize_components
from fastfish.processing.pipeline import dump,source_hash
from fastfish.processing.preview import PreviewDataset,output_grid
from fastfish.processing.stages import stage_path

class Parser(argparse.ArgumentParser):
    def error(self,message):raise ValueError(message)

def defaults(cfg):
    return dict(baseline_components=default_components(cfg),n_t_cal=cfg.get('n_t_cal',10),n_t_apply=cfg.get('n_t_apply',cfg.get('n_t',10)),baseline=cfg.get('baseline',True),protect_kms=cfg['protect_kms'],local_fit_kms=cfg['local_fit_kms'],baseline_exclusions_kms=cfg.get('baseline_exclusions_kms',[]),spectral_mode=cfg['spectral_mode'],n_b=cfg['n_b'],n_g=cfg['n_g'],native_fwhm_chan=cfg['native_fwhm_chan'],velocity_range_kms=cfg['velocity_range_kms'])

def validate_params(cfg,params):
    p=dict(params);p['baseline_components']=normalize_components(p.get('baseline_components',default_components(dict(cfg,**p))));p['n_t_cal']=smoothing_rows(p['n_t_cal']);p['n_t_apply']=smoothing_rows(p['n_t_apply'])
    if not isinstance(p['baseline'],bool):raise ValueError('baseline must be true or false')
    if not np.isfinite(p['protect_kms']) or p['protect_kms']<0:raise ValueError('Protected H I half-width must be finite and nonnegative')
    p['baseline_exclusions_kms']=normalize_windows(p['baseline_exclusions_kms'])
    if p['spectral_mode'] not in ('rebin','gaussian'):raise ValueError('spectral_mode must be rebin or gaussian')
    if int(p['n_b'])!=p['n_b'] or p['n_b']<1:raise ValueError('n_b must be a positive integer')
    if not np.isfinite([p['n_g'],p['native_fwhm_chan']]).all() or p['native_fwhm_chan']<=0 or p['n_g']<p['native_fwhm_chan']:raise ValueError('Invalid Gaussian widths')
    p['n_b']=int(p['n_b'])
    for key in ('protect_kms','local_fit_kms','n_g','native_fwhm_chan'):p[key]=float(p[key])
    p['velocity_range_kms']=list(map(float,p['velocity_range_kms']))
    output_grid(dict(cfg,**p));return p

def owned_process(pid,module,identity):
    try:
        proc=Path('/proc')/str(int(pid))/'cmdline'
        if proc.exists():
            text=proc.read_bytes().replace(b'\0',b' ').decode();return module in text and str(identity) in text
        if sys.platform!='linux':os.kill(int(pid),0);return False
    except (OSError,ValueError,TypeError):pass
    return False

class FishShell(cmd.Cmd):
    prompt='fish> '
    intro='FISH · FAST spectral reduction. Type help. Inspection is read-only; STAGE run processes the loaded selection.'
    def __init__(self,port=8767,viewer="qt"):
        super().__init__();self.viewer=viewer;self.project=None;self.path=None;self.port=port;self.ds=None;self.preview_process=None
    def emptyline(self):pass
    def default(self,line):raise ValueError('Unknown command: '+line+'. Type help.')
    def _stop_preview(self):
        if self.preview_process is not None and self.preview_process.poll() is None:self.preview_process.terminate()
        self.preview_process=None
    def onecmd(self,line):
        self.last_error=None
        try:return super().onecmd(line)
        except (ValueError,OSError,RuntimeError,KeyError,IndexError) as e:self.last_error=str(e);print('Error:',e);return False
    def _state(self):
        if self.path is None:raise ValueError('First use file in CONFIG [--selection SELECTION] [--run-root DIR]')
        return read_session(self.path)
    def _cfg(self):return json.loads(Path(self._state()['config']).read_text())
    def _data(self):
        state=self._state()
        if not state.get('selection') or not Path(state['selection']).is_file():raise ValueError('Pointing is not ready. Use go where for project imports, or prepare for configuration imports.')
        key=(state['config'],state['selection'])
        if self.ds is None or self.ds[0]!=key:self.ds=(key,PreviewDataset(*key,memory_mib=128))
        return self.ds[1]
    def _set(self,updates):
        state=self._state();params=validate_params(self._cfg(),dict(state['params'],**updates))
        write_session(self.path,dict(params=params),state['revision']);print('Settings:',json.dumps(updates))
    def do_catalog(self,arg):
        """catalog [JSON_PATH|packaged] [--profile NAME]: inspect or load an RFI search-prior catalogue; no automatic masking."""
        s=self._state()
        if not arg.strip():
            cfg=self._cfg();print('RFI catalogue:',json.dumps(s.get('rfi_catalog_info',dict(source=cfg.get('rfi_catalog','packaged'),profile=cfg.get('rfi_profile'))),indent=2));return
        p=Parser(prog='catalog',add_help=False);p.add_argument('path');p.add_argument('--profile');a=p.parse_args(shlex.split(arg))
        from fastfish.processing.session import load_session_catalog
        updated=load_session_catalog(self.path,a.path,a.profile,s['revision']);self.ds=None
        print('RFI search-prior catalogue loaded:',json.dumps(updated['rfi_catalog_info'],indent=2))
        print('Refresh the viewer. Existing products are preserved; future runs use this catalogue.')
    def do_list(self,arg):
        """list: show selected SCI metadata (UTC/MJD, timestamps, channels, frequency and size); no spectrum or pointing calculation."""
        if arg.strip():raise ValueError('Usage: list')
        state=self._state()
        if state.get('metadata'):
            records=json.loads(Path(state['metadata']).read_text())
        else:
            from fastfish.io.fits_stream import spectral_metadata
            files=state.get('input_files') or [f for i in self._cfg()['inputs'] for f in i['files']]
            records=[spectral_metadata(f) for f in files]
        from astropy.time import Time
        for number,r in enumerate(records,1):
            lo,hi=r['mjd_first_last'];utc=Time([lo,hi],format='mjd',scale='utc').iso
            first=r['freq_first_ghz']*1000;last=first+(r['channels']-1)*r['native_df_hz']/1e6
            print('{}: {}'.format(number,r['path']))
            print('  UTC ({}) {} → {} | MJD {:.9f} → {:.9f}'.format(r.get('timestamp_reference','legacy/raw'),utc[0],utc[1],lo,hi))
            print('  timestamps={}  nchan={}  frequency={:.6f}…{:.6f} MHz  df={:.9g} Hz  bytes={}'.format(r['rows'],r['channels'],first,last,r['native_df_hz'],r['bytes']))
        print('Total:',len(records),'files;',sum(r['rows'] for r in records),'timestamps. Pointing:', 'pending (go where)' if state.get('pointing_pending') else 'available' if Path(state['selection']).is_file() else 'not prepared')
    def do_project(self,arg):
        """project DIRECTORY: establish org_spec, org_cal and feed_info locations without reading data."""
        words=shlex.split(arg)
        if len(words)!=1:raise ValueError('Usage: project DIRECTORY')
        root=Path(words[0]).expanduser().resolve()
        if not (root/'org_spec').is_dir():raise ValueError('Missing science directory: '+str(root/'org_spec'))
        self._stop_preview();self.project=root;self.path=None;self.ds=None
        print('Project:',root,'— open selected FITS files, then go where')
    def _project_files(self,a):
        import tempfile
        from fastfish.io.file_input import expand_inputs
        from fastfish.io.project_input import register
        files=expand_inputs(a.inputs,self.project);cfg,records=register(self.project,files)
        if a.run_root:
            root=Path(a.run_root).expanduser().resolve()
            if root.exists() and any(root.iterdir()):raise ValueError('Run directory is not empty')
            root.mkdir(parents=True,exist_ok=True)
        else:
            parent=self.project/'fish_runs'/'file_subsets';parent.mkdir(parents=True,exist_ok=True);root=Path(tempfile.mkdtemp(prefix='files-',dir=str(parent))).resolve()
        cfg['output_dir']=str(root/'coordinate_setup');config=root/'config.json';dump(config,cfg);dump(root/'science_metadata.json',records)
        path=root/'fish_session.json';write_session(path,dict(config=str(config),selection=str(root/'selection.json'),root=str(root),params=validate_params(cfg,defaults(cfg)),index=0,input_files=[str(f) for f in files],project_root=str(self.project),metadata=str(root/'science_metadata.json'),pointing_pending=True))
        self._stop_preview();self.path=path;self.ds=None
        print('Opened',len(files),'SCI files;',sum(r['rows'] for r in records),'timestamps. Pointing not calculated. Use go where. Session:',path)
    def do_file(self,arg):
        """file in FITS [FITS ...] | file in /list LISTNAME [--config CONFIG]. Also accepts an observation configuration/example."""
        from fastfish.io.file_input import native_name,expand_inputs,subset_config,cached_subset
        import tempfile
        p=Parser(prog='file in',add_help=False);p.add_argument('direction',choices=['in']);p.add_argument('inputs',nargs='+');p.add_argument('--config');p.add_argument('--selection');p.add_argument('--run-root');a=p.parse_args(shlex.split(arg))
        direct=a.inputs[0]=='/list' or any(native_name(x) for x in a.inputs)
        if not direct:
            if len(a.inputs)!=1 or a.config:raise ValueError('Use one observation config, or FITS files with --config CONFIG')
            return self._open_config(arg)
        if self.project is not None and not a.config:
            if a.selection:raise ValueError('Project imports compute their full footprint; --selection is not used')
            return self._project_files(a)
        files=expand_inputs(a.inputs);source=None;seed=None;inherited=None
        if a.config:
            source=Path(a.config).expanduser()
            if not source.is_file():source=PROJECT_ROOT/'Example'/a.config/'config.json'
        elif self.path:
            state=self._state();source=Path(state.get('observation_config') or state['config']);seed=Path(state.get('observation_selection') or state['selection']);inherited=state['params']
        else:
            candidates=[]
            for candidate in (PROJECT_ROOT/'Example').glob('*/config.json'):
                try:subset_config(candidate,files);candidates.append(candidate)
                except (ValueError,OSError,KeyError):pass
            if len(candidates)!=1:raise ValueError('Cannot uniquely associate these files with calibration and pointing. Use --config CONFIG or first file in OBSERVATION.')
            source=candidates[0]
        cfg=subset_config(source,files)
        if seed is None:
            original=json.loads(source.read_text());value=a.selection or original.get('interactive_selection')
            if value:seed=(source.parent/Path(value).expanduser()).resolve()
        selected=cached_subset(seed,files) if seed and seed.is_file() else None
        # A new session prevents a small-file test from overwriting a full observation run.
        if a.run_root:
            root=Path(a.run_root).expanduser().resolve()
            if root.exists() and any(root.iterdir()):raise ValueError('Subset run directory is not empty; choose a new --run-root')
            root.mkdir(parents=True,exist_ok=True)
        else:
            parent=Path(self._state()['root']).parent/'file_subsets' if self.path else Path.cwd()/'fish_runs'/'file_subsets'
            parent.mkdir(parents=True,exist_ok=True);root=Path(tempfile.mkdtemp(prefix='files-',dir=str(parent))).resolve()
        cfg['output_dir']=str(root/'coordinate_setup');cfg['fish_run_root']=str(root);cfg.pop('interactive_selection',None);config=root/'config.json';dump(config,cfg)
        params=validate_params(cfg,inherited or defaults(cfg));path=root/'fish_session.json'
        if selected is not None:dump(root/'selection.json',selected)
        write_session(path,dict(config=str(config),selection=str(root/'selection.json'),root=str(root),params=params,index=0,input_files=[str(f) for f in files],observation_config=str(source.resolve()),observation_selection=str(seed) if seed else None))
        self._stop_preview();self.path=path;self.ds=None
        if selected is None:self.do_prepare('')
        d=self._data();index=d.closest(cfg['header']['CRVAL1'],cfg['header']['CRVAL2']);write_session(path,dict(index=index))
        print('Imported',len(files),'FITS files;',len(d.ra),'selected spectra. Session:',path)
    def _open_config(self,arg):
        """file in CONFIG_OR_EXAMPLE [--selection FILE] [--run-root DIR]: open an observation, not an individual uncalibrated FITS file."""
        p=Parser(prog='file in',add_help=False);p.add_argument('direction',choices=['in']);p.add_argument('config');p.add_argument('--selection');p.add_argument('--run-root');a=p.parse_args(shlex.split(arg))
        config=Path(a.config).expanduser()
        if not config.is_file():config=PROJECT_ROOT/'Example'/a.config/'config.json'
        config=config.resolve();cfg=json.loads(config.read_text());root=Path(a.run_root or cfg.get('fish_run_root',str(Path.cwd()/'fish_runs'/cfg.get('name',config.parent.name)))).expanduser().resolve()
        selected=a.selection or cfg.get('interactive_selection') or str(config.parent/cfg.get('output_dir','output')/'selection.json')
        selection=Path(selected).expanduser();selection=(config.parent/selection).resolve() if not selection.is_absolute() else selection.resolve()
        new_path=root/'fish_session.json'
        self.project=None
        write_session(new_path,dict(config=str(config),selection=str(selection),root=str(root),params=validate_params(cfg,defaults(cfg)),index=0))
        self._stop_preview();self.path=new_path;self.ds=None
        if selection.is_file():
            d=self._data();index=d.closest(cfg['header']['CRVAL1'],cfg['header']['CRVAL2']);write_session(self.path,dict(index=index));print('Loaded',cfg.get('name'),len(d.ra),'spectra. Session:',self.path)
        else:print('Configuration loaded. Selection missing; use prepare. Session:',self.path)
    def do_prepare(self,arg):
        """prepare: calculate inventory/pointing and the spatial selection in this session's own directory."""
        if arg.strip():raise ValueError('Usage: prepare')
        if self._state().get('project_root'):raise ValueError('Use go where to calculate project-file pointing')
        s=self._state();cfg=self._cfg();base=Path(s['config']).parent;root=Path(s['root'])
        for key in ('input_root','extracted_cal_root','tcal_path'):
            if cfg.get(key):cfg[key]=str((base/Path(cfg[key]).expanduser()).resolve())
        if cfg.get('rfi_catalog','packaged')!='packaged':cfg['rfi_catalog']=str((base/cfg['rfi_catalog']).resolve())
        if cfg.get('pointing'):cfg['pointing']={k:str((base/v).resolve()) for k,v in cfg['pointing'].items()}
        for item in cfg.get('inputs',[]):
            item['reference']=str((base/item['reference']).resolve());item['files']=[str((base/f).resolve()) for f in item['files']]
        cfg['output_dir']=str(root/'coordinate_setup');snapshot=root/'input_config.json';dump(snapshot,cfg)
        from fastfish.commands.prepare_trial import prepare
        prepare(snapshot,root);write_session(self.path,dict(config=str(snapshot),selection=str(root/'selection.json')));self.ds=None;print('Prepared',len(self._data().ra),'spectra')
    def do_go(self,arg):
        """go where [/workers 1..64]: calculate and display selected-file pointing; default 2 workers."""
        words=shlex.split(arg)
        if words!=['where'] and not (len(words)==3 and words[:2]==['where','/workers']):raise ValueError('Usage: go where [/workers 1..64]')
        workers=int(words[2]) if len(words)==3 else None
        if workers is not None and not 1<=workers<=64:raise ValueError('Pointing workers must be from 1 to 64')
        state=self._state()
        if state.get('pointing_pending'):
            from fastfish.io.project_input import pointing
            cfg=pointing(state['config'],state['metadata'],state['root'],workers=workers);write_session(self.path,dict(pointing_pending=False),state['revision']);self.ds=None
            print('Footprint:',len(self._data().ra),'supported timestamps; no spatial or calibration crop')
        self._view('where')
    def do_get(self,arg):
        """get INDEX | next | previous: choose a spectrum for plots."""
        d=self._data();s=self._state();value=arg.strip()
        index=s['index']+1 if value=='next' else s['index']-1 if value=='previous' else int(value)
        if not 0<=index<len(d.ra):raise ValueError('Spectrum index outside selection')
        write_session(self.path,dict(index=index));print('Spectrum',index,'RA/Dec',float(d.ra[index]),float(d.dec[index]))
    def do_set(self,arg):
        """set KEY JSON_VALUE: change settings, e.g. set n_b 4; set velocity_range_kms [-100,100]; set baseline_exclusions_kms [[700,900]]."""
        key,sep,value=arg.strip().partition(' ')
        if not sep or key not in self._state()['params']:raise ValueError('Use set KEY JSON_VALUE; settings lists the keys')
        try:parsed=json.loads(value)
        except json.JSONDecodeError:
            if key=='spectral_mode':parsed=value
            else:raise ValueError('Value must be JSON, e.g. 10, true, [-100,100]')
        self._set({key:parsed})
    def do_settings(self,arg):
        """settings: show the settings currently shared with the plot window."""
        print(json.dumps(self._state(),indent=2))
    def do_tsyscal(self,arg):
        """tsyscal [N_CAL N_APPLY] | tsyscal run: inspect Tsys and its two smoothing windows, or generate references for all loaded spectra."""
        words=shlex.split(arg)
        if words==['run']:return self._run('reference')
        if words:
            if len(words)!=2:raise ValueError('Usage: tsyscal [N_CAL N_APPLY] or tsyscal run')
            self._set(dict(n_t_cal=int(words[0]),n_t_apply=int(words[1])))
        self._view('reference')
    def do_applycal(self,arg):
        """applycal [run]: inspect calibrated AA/BB, or generate gain-transfer records for all loaded spectra."""
        self._stage('calibration',arg)
    def do_baseline(self,arg):
        """baseline [run]: inspect/edit independent AA/BB fits, or fit all loaded spectra with the current windows."""
        self._stage('baseline',arg)
    def do_table(self,arg):
        """table [run]: inspect resampling, or build the disk-backed final spectral table and table_manifest.json for gridding."""
        self._stage('resample',arg)
    def do_xy_map(self,arg):
        """xy_map run [--validate]: grid and combine the prepared table. xy_map alone reports the current cube location."""
        words=shlex.split(arg)
        if not words:
            s=self._state();path=stage_path(s['root'],'cube',self._data(),s['params']);print('Cube directory:',path,'(ready)' if (path/'summary.json').is_file() else '(not produced; use xy_map run)');return
        if words not in (['run'],['run','--validate']):raise ValueError('Usage: xy_map run [--validate]')
        self._run('cube',validate='--validate' in words)
    def _stage(self,stage,arg):
        if arg.strip()=='run':self._run(stage)
        elif not arg.strip():self._view(stage)
        else:raise ValueError('Use the command alone to inspect, or add run to process the loaded selection')
    def _view(self,view):
        if self.viewer=='none':return
        s=self._state();self._data()
        if self.viewer=='qt':
            if os.environ.get('SSH_CONNECTION') or (sys.platform=='linux' and not os.environ.get('DISPLAY')):
                print('Open the native viewer on your computer: fish-viewer --remote USER@HOST --session-state '+shlex.quote(str(self.path))+(' --view where' if view=='where' else ''))
                return
            if self.preview_process is None or self.preview_process.poll() is not None:
                self.preview_process=subprocess.Popen([sys.executable,'-m','fastfish.qt_viewer','--session-state',str(self.path),'--view',view])
            return
        if self.preview_process is None or self.preview_process.poll() is not None:
            with socket.socket() as test:
                try:test.bind(('127.0.0.1',self.port))
                except OSError:raise ValueError('Preview port is occupied; restart fish with --port PORT. No existing service was changed.')
            root=Path(s['root']);log=root/'fish_preview.log'
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
            with log.open('ab') as output:self.preview_process=subprocess.Popen([sys.executable,'-m','fastfish.interactive_app','--session-state',str(self.path),'--port',str(self.port)],stdout=output,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            dump(root/'fish_preview.pid',dict(pid=self.preview_process.pid,session=str(self.path)))
        url='http://127.0.0.1:{}/?view={}&revision={}'.format(self.port,view,s['revision']);print('Plot:',url)
        print('For an SSH session, forward this same port to your computer. Plot edits sync to this prompt; reopening a plot refreshes terminal edits.')
        if not os.environ.get('SSH_CONNECTION') and (sys.platform=='darwin' or os.environ.get('DISPLAY')):webbrowser.open(url)
    def _run(self,stage,validate=False):
        s=self._state();d=self._data();params=validate_params(d.cfg,s['params']);root=Path(s['root'])
        dependency={'calibration':'reference','baseline':'calibration','resample':'baseline','cube':'resample'}.get(stage)
        if dependency and not (stage_path(root,dependency,d,params)/'summary.json').is_file():raise ValueError('Complete '+{'reference':'tsyscal','calibration':'applycal','baseline':'baseline','resample':'table'}[dependency]+' run first with these settings')
        with (root/'.launch.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            job=root/'active_job.json'
            if job.exists() and owned_process(json.loads(job.read_text()).get('pid'),'fastfish.native_stages',root/'jobs'):raise ValueError('A job is already running; use status or wait')
            jobs=root/'jobs';jobs.mkdir(exist_ok=True);key=str(time.time_ns());recipe=jobs/(key+'.json');log=jobs/(key+'.log')
            dump(recipe,dict(config=s['config'],selection=s['selection'],root=str(root),stage=stage,params=params,source_sha256=source_hash(),catalog_sha256=d.catalog_sha256,validate=validate))
            with log.open('ab') as output:p=subprocess.Popen([sys.executable,'-m','fastfish.native_stages','--recipe',str(recipe)],stdout=output,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
            dump(job,dict(pid=p.pid,stage=stage,recipe=str(recipe),log=str(log)))
        print('Started',stage,'for',len(d.ra),'spectra; PID',p.pid,'log',log)
    def _finish_table(self):
        s=self._state();d=self._data();folder=stage_path(s['root'],'resample',d,s['params'])
        if (folder/'summary.json').is_file():
            dump(Path(s['root'])/'table_manifest.json',dict(format='fastfish-gridding-table-v1',source_sha256=source_hash(),stage=str(folder),params=s['params'],selection=s['selection'],frequency_ghz=output_grid(dict(d.cfg,**s['params'])).tolist(),processed=str(folder/'processed.json'),storage='Disk-backed arrays plus RA/Dec/MJD/Doppler metadata; no native spectra copied'))
    def do_status(self,arg):
        """status: show running job, stage checkpoints and the table/cube locations."""
        s=self._state();root=Path(s['root']);d=self._data()
        for stage in ('reference','calibration','baseline','resample','cube'):
            p=stage_path(root,stage,d,s['params']);print(stage, 'complete' if (p/'summary.json').is_file() else 'pending',p)
        job=root/'active_job.json'
        if job.is_file():
            record=json.loads(job.read_text());print('Job:',record['stage'],'running' if owned_process(record['pid'],'fastfish.native_stages',root/'jobs') else 'not running','log:',record['log'])
        if (root/'status.json').is_file():print('Progress:',(root/'status.json').read_text().strip())
        self._finish_table()
    def do_wait(self,arg):
        """wait [SECONDS]: wait for the current background stage; Ctrl-C returns to the prompt without stopping it."""
        s=self._state();root=Path(s['root']);path=root/'active_job.json'
        if not path.is_file():print('No job');return
        job=json.loads(path.read_text());deadline=time.monotonic()+float(arg) if arg.strip() else None
        while owned_process(job['pid'],'fastfish.native_stages',root/'jobs'):
            if deadline and time.monotonic()>deadline:print('Still running');return
            try:time.sleep(.25)
            except KeyboardInterrupt:print('Job continues; use stop to terminate it.');return
        self._finish_table();recipe=json.loads(Path(job['recipe']).read_text());p=stage_path(root,job['stage'],self._data(),recipe['params'])
        if (p/'summary.json').is_file():print('Complete:',job['stage'],p)
        else:raise RuntimeError('Stage stopped or failed. Inspect '+job['log'])
    def do_stop(self,arg):
        """stop: terminate only the identified active science job in this run directory."""
        s=self._state();p=Path(s['root'])/'active_job.json'
        if p.is_file():
            job=json.loads(p.read_text())
            if owned_process(job['pid'],'fastfish.native_stages',Path(s['root'])/'jobs'):os.killpg(job['pid'],signal.SIGTERM);print('Stop requested');return
        print('No active owned job')
    def do_save(self,arg):
        """save FILE: save observation paths, settings and spectrum index as a small JSON session."""
        if not arg.strip():raise ValueError('Usage: save FILE')
        dump(Path(shlex.split(arg)[0]).expanduser(),self._state());print('Session saved')
    def do_restore(self,arg):
        """restore FILE: restore a saved JSON session, keeping its original run directory."""
        if len(shlex.split(arg))!=1:raise ValueError('Usage: restore FILE')
        state=read_session(Path(shlex.split(arg)[0]).expanduser());cfg=json.loads(Path(state['config']).read_text());state['params']=validate_params(cfg,state['params'])
        new_path=Path(state['root'])/'fish_session.json';write_session(new_path,{k:v for k,v in state.items() if k not in ('revision','format')});self._stop_preview();self.path=new_path;self.ds=None;self.project=Path(state['project_root']) if state.get('project_root') else None;print('Session restored:',self.path)
    def do_exit(self,arg):
        """exit: close this shell's plot service; background science jobs continue."""
        self._stop_preview()
        return True
    do_quit=do_exit
    def do_EOF(self,arg):print();return self.do_exit(arg)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--viewer',choices=['browser','qt','none'],default='qt',help='Plot viewer (default: qt); none suppresses windows for unattended runs');p.add_argument('--port',type=int,default=8767);p.add_argument('--execute',help='Read commands from a text file, useful for reproducible reductions');a=p.parse_args()
    if not 1<=a.port<=65535:p.error('port must be between 1 and 65535')
    shell=FishShell(a.port,a.viewer)
    if a.execute:
        for line in Path(a.execute).read_text().splitlines():
            line=line.strip()
            if line and not line.startswith('#'):
                print('fish>',line)
                if shell.onecmd(line):break
                if shell.last_error:shell.do_exit('');return 1
        shell.do_exit('')
    else:
        while True:
            try:shell.cmdloop();break
            except KeyboardInterrupt:print('\nUse exit to leave, or stop to stop a science job.')

if __name__=='__main__':raise SystemExit(main())
