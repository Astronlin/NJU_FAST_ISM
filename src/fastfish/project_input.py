"""Lazy project discovery: SCI metadata first, telemetry only for go where."""
from pathlib import Path
import json,re
import numpy as np
from .fits_stream import spectral_metadata
from .pipeline import dump,_auto_header


def processing_defaults():
    return dict(rest_freq_ghz=1.420405751768,native_df_hz=476.837158203125,
        coordinate_frame='icrs',spectral_mode='rebin',n_b=2,n_g=2.,native_fwhm_chan=1.,baseline=True,
        local_fit_kms=1500.,protect_kms=250.,alternative_fit_kms=1000.,gain_degree=3,
        chunk_rows=32,grid_channel_block=64,grid_row_block=256,workers=64,io_workers=4,
        velocity_range_kms=[-250.,250.],pixel_arcmin=1.,spatial_margin_pixels=1.910974050648043,
        n_t_apply=10,n_t_cal=10,science_timestamp_reference='start',baseline_exclusions_kms=[],rfi_catalog='packaged',rfi_profile=None)


def register(project,files):
    project=Path(project).resolve();records=[spectral_metadata(p) for p in files]
    first=records[0]
    for r in records:
        if r['channels']!=first['channels'] or not np.isclose(r['native_df_hz'],first['native_df_hz'],rtol=1e-8):raise ValueError('Selected files have incompatible native spectral grids')
    cfg=processing_defaults();cfg.update(name=project.name,project_root=str(project),input_root=str(project/'org_spec'),
        extracted_cal_root=str(project/'org_cal'),native_df_hz=first['native_df_hz'],inputs=[],scans=[],pointing={},header=None)
    candidates=list((Path(__file__).resolve().parents[2]/'Tcal').glob('*.fits'))
    cfg['tcal_path']=str(candidates[0]) if len(candidates)==1 else None
    return cfg,records


def associations(project,records):
    project=Path(project);sheets=sorted((project/'feed_info').glob('*.xlsx'));cal=project/'org_cal'
    # File/directory names identify candidates; telemetry content is read only in pointing().
    scans={p.name for p in cal.iterdir() if p.is_dir()} if cal.is_dir() else set()
    scans.update(re.sub(r'_\d{4}_\d{2}_\d{2}.*$','',p.stem) for p in sheets)
    items={};mapping={}
    for r in records:
        name=Path(r['path']).name;beam=re.search(r'-(M\d{2})_',name)
        choices=[s for s in scans if name.startswith(s+'_') or name.startswith(s+'-')]
        if not beam or len(choices)!=1:raise ValueError('Cannot uniquely identify observation/beam for '+r['path'])
        scan=choices[0];b=beam.group(1);matches=[p for p in sheets if p.stem==scan or p.stem.startswith(scan+'_')]
        if len(matches)!=1:raise ValueError('Expected one telemetry sheet for '+scan+'; candidates: '+', '.join(map(str,matches)))
        mapping[scan]=str(matches[0]);key=(scan,b)
        if key not in items:items[key]=dict(scan=scan,beam=b,reference=str(project/'org_cal'/scan/b/'Tsys_reference_not_precomputed.fits'),files=[])
        items[key]['files'].append(r['path'])
    return list(items.values()),mapping


# Workers receive paths, never SCI spectral arrays or large pickled telemetry tables.
_POINTING_TABLES={}

def _pointing_init():
    global _POINTING_TABLES
    _POINTING_TABLES={}
    from .time_tables import initialize_earth_orientation
    initialize_earth_orientation()

def _pointing_file(task):
    from .xyz2celestial import reconstruct_pointing
    import pandas as pd
    scan,beam,path,table_path,out,cache,time_reference=task
    if scan not in _POINTING_TABLES:
        with np.load(table_path) as z:_POINTING_TABLES[scan]=pd.DataFrame({k:z[k] for k in z.files})
    entry,t=spectral_metadata(path,read_times=True,time_reference=time_reference)
    ra,dec=reconstruct_pointing(t,_POINTING_TABLES[scan],beam=int(beam[1:]));ids=np.flatnonzero(np.isfinite(ra)&np.isfinite(dec))
    if len(ids):
        # Publish directly to the one project cache; no second session coordinate copy.
        import tempfile,os
        from .coordinate_cache import publish
        folder=Path(cache['path']).parent;folder.mkdir(parents=True,exist_ok=True)
        fd,temp=tempfile.mkstemp(prefix='.pending-',suffix='.npz',dir=str(folder));os.close(fd)
        try:
            np.savez_compressed(temp,ids=ids,mjd=t[ids],ra=ra[ids],dec=dec[ids],doppler_pending=np.array(True),time_reference=np.array('midpoint'),raw_time_reference=np.array(time_reference))
            publish(cache,temp,len(ids),len(t)-len(ids))
        finally:
            if os.path.exists(temp):os.unlink(temp)
        out=cache['path']
    import resource
    return dict(entry=entry,coordinates=str(out),kept=len(ids),excluded=len(t)-len(ids),peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _ordered_results(pool,tasks,workers):
    """Limit queued files and cancel work not started on failure/interruption."""
    from collections import deque
    pending=deque();source=iter(tasks)
    for _ in range(2*workers):
        task=next(source,None)
        if task is not None:pending.append(pool.submit(_pointing_file,task))
    try:
        while pending:
            future=pending.popleft();yield future.result()
            task=next(source,None)
            if task is not None:pending.append(pool.submit(_pointing_file,task))
    finally:
        for future in pending:future.cancel()

def pointing(config_path,metadata_path,root,workers=None):
    import fcntl
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    with (root/'.pointing.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Pointing is already running in this session')
        return _pointing_unlocked(config_path,metadata_path,root,workers)


def _pointing_unlocked(config_path,metadata_path,root,workers=None):
    import time
    started=time.perf_counter()
    from astropy.table import Table
    from astropy.utils import iers
    import pandas as pd
    from .telemetry_fast import telemetry_cache
    from .xyz2celestial import reconstruct_pointing
    from .interactive_core import output_grid
    from .rfi_catalog import load_catalog
    from astropy.time import Time
    cfg=json.loads(Path(config_path).read_text());records=json.loads(Path(metadata_path).read_text());root=Path(root)
    reference=cfg.get('science_timestamp_reference','start')
    if any(r.get('raw_timestamp_reference')!=reference for r in records):
        records=[spectral_metadata(r['path'],time_reference=reference) for r in records];dump(metadata_path,records)
    items,sheets=associations(cfg['project_root'],records);scans=sorted(sheets);cfg.update(inputs=items,scans=scans,pointing=sheets)
    iers.conf.auto_download=False;iers.conf.auto_max_age=None
    iers.earth_orientation_table.set(iers.IERS_B(Table.read(Path(__file__).with_name('iers_b.fits'))))
    import os,time,multiprocessing
    from concurrent.futures import ProcessPoolExecutor
    workers=cfg.get('pointing_workers',2) if workers is None else workers
    if isinstance(workers,bool) or int(workers)!=workers or not 1<=workers<=64:raise ValueError('Pointing workers must be an integer from 1 to 64')
    requested=int(workers)
    from .coordinate_cache import context,candidate,lookup
    common=context(cfg['project_root'],sheets);by_path={r['path']:r for r in records};cached={};cache_info={};ordered=[]
    for item in items:
        for path in item['files']:
            index=len(ordered);ordered.append((item,path));cache=candidate(cfg['project_root'],item['scan'],item['beam'],by_path[path],common);cache_info[index]=cache;hit=lookup(cache)
            if hit is not None:cached[index]=dict(hit,entry=by_path[path])
    missing=[i for i in range(len(ordered)) if i not in cached];needed_scans={ordered[i][0]['scan'] for i in missing};tables={}
    for scan in sorted(needed_scans):
        entries=[by_path[ordered[i][1]] for i in missing if ordered[i][0]['scan']==scan];bounds=[v for r in entries for v in r['mjd_first_last']]
        info=telemetry_cache(scan,sheets[scan],min(bounds),max(bounds),root/'pointing_cache');tables[scan]=info['path']
    telemetry_seconds=time.perf_counter()-started;tasks=[]
    for i in missing:
        item,path=ordered[i];out=root/'coordinates'/item['scan']/item['beam']/(Path(path).stem+'.npz')
        tasks.append((item['scan'],item['beam'],path,tables[item['scan']],str(out),cache_info[i],cfg.get('science_timestamp_reference','start')))
    workers=min(requested,len(tasks)) if tasks else 0
    selected=[];all_ra=[];all_dec=[];excluded=0;computed=[]
    def collect(iterator):
        for number,result in enumerate(iterator,1):
            computed.append(result)
            print('POINTING',number,'/',len(tasks),Path(result['entry']['path']).name,result['kept'],'timestamps',flush=True)
    print('COORDINATE CACHE',len(cached),'reused;',len(tasks),'files to calculate',flush=True)
    if workers==1:
        _pointing_init();collect(map(_pointing_file,tasks))
    elif workers>1:
        names=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS');previous={k:os.environ.get(k) for k in names}
        try:
            for k in names:os.environ[k]='1'
            with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),initializer=_pointing_init) as pool:
                collect(_ordered_results(pool,tasks,workers))
        finally:
            for k,v in previous.items():
                if v is None:os.environ.pop(k,None)
                else:os.environ[k]=v
    by_index=dict(cached);by_index.update(zip(missing,computed));results=[by_index[i] for i in range(len(ordered))]
    position=0
    for item in items:
        contexts=[]
        for path in item['files']:
            result=results[position];position+=1;excluded+=result['excluded']
            if not result['kept']:continue
            contexts.append(dict(entry=result['entry'],coordinates=result['coordinates']))
            with np.load(result['coordinates']) as z:all_ra.append(z['ra']);all_dec.append(z['dec'])
        if contexts:selected.append(dict(item,files=contexts))
    if not selected:raise ValueError('No selected timestamps are covered by valid telemetry')
    cfg['header']=cfg.get('header') or _auto_header(np.concatenate(all_ra),np.concatenate(all_dec),cfg['pixel_arcmin'],cfg.get('coordinate_frame','icrs'))
    catalogue=load_catalog('packaged',Path(config_path).parent);dates=[Time(v,format='mjd').iso[:10] for r in records for v in r['mjd_first_last']];profiles=[]
    for name,p in catalogue['profiles'].items():
        bounds=p.get('observation_date_range',[])
        if set(scans)<=set(p['scan_priors']) and len(bounds)==2 and bounds[0]<=min(dates)<=max(dates)<=bounds[1]:profiles.append(name)
    cfg['rfi_profile']=profiles[0] if len(profiles)==1 else None
    # Commit the resolved configuration only after the full pointing calculation succeeds.
    dump(config_path,cfg);dump(root/'selection.json',dict(items=selected,header=cfg['header'],output_freq_ghz=output_grid(cfg).tolist(),spectra=sum(len(x) for x in all_ra),excluded_rows=excluded,spatial_crop=False,calibration_time_crop=False))
    dump(root/'pointing_performance.json',dict(workers=workers,files=len(records),cache_hits=len(cached),computed_files=len(tasks),telemetry_seconds=telemetry_seconds,total_seconds=time.perf_counter()-started,max_worker_peak_rss_kib=max(r['peak_rss_kib'] for r in results),memory_note='Maximum per-worker high-water RSS, not simultaneous total RAM; Linux KiB'))
    import shutil
    shutil.rmtree(str(root/'pointing_cache'),ignore_errors=True)
    return cfg
