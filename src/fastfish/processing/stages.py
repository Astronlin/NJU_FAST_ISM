"""Compact restartable batch records for the four interactive processing stages."""
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import multiprocessing
import time
import traceback
from fastfish.processing.scheduler import memory_map,worker_estimate
import numpy as np
from scipy.interpolate import interp1d
from fastfish.processing.preview import (PreviewDataset,calibrate,reconstruct_calibrated,
    fit_polarization_baselines,reconstruct_polarization_baselines,resample,output_grid)
from fastfish.calibration.reference import reference_data
from fastfish.calibration.extracted import read_reference
from fastfish.io.fits_stream import fixed_table,contiguous_rows
from fastfish.processing.pipeline import dump,read,source_hash,grid,combine,validate
from fastfish.rfi.detection import C
from fastfish.baseline.components import fit_components,reconstruct_components
from fastfish.rfi.catalog import load_catalog

_WORKER=None


def fingerprint(stage,ds,params):
    record=dict(input=ds.input_signature,source=source_hash(),catalog=ds.catalog_sha256,stage=stage,n_t_apply=params.get('n_t_apply',10),n_t_cal=params.get('n_t_cal',10))
    if stage in ('baseline','resample','cube'):
        record['baseline']={k:params[k] for k in ('baseline','protect_kms','local_fit_kms')}
        record['baseline']['baseline_exclusions_kms']=params.get('baseline_exclusions_kms',[])
        record['baseline']['baseline_components']=params.get('baseline_components')
    if stage in ('resample','cube'):
        record['spectral']={k:params[k] for k in ('spectral_mode','n_b','n_g','native_fwhm_chan','velocity_range_kms')}
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()


def stage_path(root,stage,ds,params):
    return Path(root)/'stages'/stage/fingerprint(stage,ds,params)[:16]


def init_worker(config,selection,params,root,stage):
    global _WORKER
    ds=PreviewDataset(config,selection,memory_mib=128)
    _WORKER=(ds,params,Path(root),stage)


def _groups(ds):
    groups={}
    for fi,ctx in enumerate(ds.files):groups.setdefault((ctx['scan'],ctx['beam']),[]).append(fi)
    return groups


def run_group(key):
    ds,params,root,stage=_WORKER;cfg=dict(ds.cfg,**params);group=_groups(ds)[tuple(key)]
    out=stage_path(root,stage,ds,params)/key[0]/key[1];out.mkdir(parents=True,exist_ok=True)
    done=out/'done.json'
    if done.exists():
        record=read(done)
        if all(Path(p).is_file() for p in record['outputs']):return record
    calfolder=stage_path(root,'calibration',ds,params)/key[0]/key[1]
    basefolder=stage_path(root,'baseline',ds,params)/key[0]/key[1]
    reffolder=stage_path(root,'reference',ds,params)/key[0]/key[1]
    products=[];outputs=[];rows=0
    if stage=='reference':
        ctx=ds.files[group[0]];f=ds.frequency(group[0]);times=[t for fi in group for t in (ds.files[fi]['coords']['mjd'][0],ds.files[fi]['coords']['mjd'][-1])]
        pack=read_reference(ctx,min(times),max(times),f,cfg)
        p=out/'temperature_reference.npz';np.savez_compressed(p,mjd=pack['mjd'],frequency_ghz=pack['frequency_ghz'],smoothed=pack['smoothed'])
        dump(out/'reference_provenance.json',pack['info'])
        outputs=[str(p),str(out/'reference_provenance.json')];rows=sum(len(ds.files[fi]['coords']['ids']) for fi in group)
    else:
        if stage=='calibration':
            with np.load(reffolder/'temperature_reference.npz') as p:reference={k:p[k] for k in p.files}
            ti=interp1d(reference['mjd'],reference['smoothed'],axis=0,bounds_error=True)
        for fi in group:
            ctx=ds.files[fi];c=ctx['coords'];f=ds.frequency(fi);n=len(c['ids']);stem=Path(ctx['entry']['path']).stem
            nchunks=int(np.ceil(n/cfg['chunk_rows']));chunk_records=[]
            if stage in ('baseline','resample'):
                with np.load(calfolder/(stem+'_calibration.npz')) as p:
                    coefficient=p['coefficients'];edges=p['reference_edges_ghz'];factors=p['frequency_scale'].copy() if 'frequency_scale' in p.files else ds.scales(fi)
                masks=read(calfolder/(stem+'_intervals.json'))
            if stage=='resample':
                with np.load(basefolder/(stem+'_baseline.npz')) as p:base={k:p[k] for k in p.files}
                q=output_grid(cfg)
                fields=('data','delta','rfi','sensitivity_missing') if cfg.get('baseline_components') is not None else ('data','delta','rfi')
                product={k:str(out/(stem+'_'+k+'.npy')) for k in fields};product['coordinates']=ctx['coordinates']
                buffers={k:np.lib.format.open_memmap(product[k],mode='w+',dtype='float32',shape=(n,len(q))) for k in fields}
            if stage=='calibration':factors=ds.scales(fi)
            coefficients=[];global_co=[];local_co=[];alt_co=[];periods=[];reports=[];component_chunks=[]
            streaming=contiguous_rows(c['ids'])
            with fixed_table(ctx['entry']['path'],random_access=not streaming,sequential_access=streaming) as (tab,_):
                for chunk in range(nchunks):
                    a=chunk*cfg['chunk_rows'];b=min(n,a+cfg['chunk_rows']);sl=slice(a,b);fac=factors[sl]
                    raw=tab['DATA'][c['ids'][sl],:,:2].astype(float)
                    if stage=='calibration':
                        result=calibrate(raw,f,ti(c['mjd'][sl]),reference['frequency_ghz'],fac,ds.priors[ctx['scan']],cfg,compact=True)
                        coefficients.append(result['coefficients']);edges=result['reference_edges_ghz']
                        chunk_records.append(dict(start=a,stop=b,intervals=result['intervals'],tracks=result['tracks']))
                    elif stage=='baseline':
                        pol=reconstruct_calibrated(raw,f,coefficient[sl],edges)
                        if cfg.get('baseline_components') is not None:
                            components,report=fit_components(pol,f,fac,masks[chunk]['intervals'],cfg);component_chunks.append(json.dumps(components));reports.append(report);continue
                        fit=fit_polarization_baselines(pol,f,fac,masks[chunk]['intervals'],cfg,sensitivity=True,compact=True)
                        global_co.append(fit['global_coeff']);local_co.append(fit['local_coeff']);alt_co.append(fit['alternative_coeff'])
                        periods.append([fit['report']['polarizations'][p]['period_MHz'] for p in ('AA','BB')]);reports.append(fit['report'])
                    else:
                        pol=reconstruct_calibrated(raw,f,coefficient[sl],edges)
                        if 'component_records_json' in base:
                            components=json.loads(str(base['component_records_json'][chunk]));residual_pol,_=reconstruct_components(pol,f,fac,cfg,components);residual=residual_pol.mean(axis=2);delta_native=np.zeros_like(residual)
                            buffers['sensitivity_missing'][sl]=1.
                        else:
                            residual,delta_native=reconstruct_polarization_baselines(pol,f,fac,base['global_coeff'][sl],base['local_coeff'][sl],base['period_MHz'][chunk],cfg,base['alternative_coeff'][sl])
                        values,fraction=resample(residual,f,fac,masks[chunk]['intervals'],cfg,q)
                        buffers['data'][sl]=values;buffers['rfi'][sl]=fraction
                        buffers['delta'][sl]=resample(delta_native,f,fac,masks[chunk]['intervals'],cfg,q)[0]
            if stage=='calibration':
                p=out/(stem+'_calibration.npz');np.savez_compressed(p,coefficients=np.concatenate(coefficients),reference_edges_ghz=edges,frequency_scale=factors)
                masks_path=out/(stem+'_intervals.json');dump(masks_path,chunk_records);outputs.extend([str(p),str(masks_path)])
            elif stage=='baseline':
                p=out/(stem+'_baseline.npz')
                if cfg.get('baseline_components') is not None:np.savez_compressed(p,component_records_json=np.asarray(component_chunks),baseline_mode='ordered_components')
                else:np.savez_compressed(p,global_coeff=np.concatenate(global_co),local_coeff=np.concatenate(local_co),alternative_coeff=np.concatenate(alt_co),period_MHz=periods,baseline_mode='per_polarization')
                report_path=out/(stem+'_fit_reports.json');dump(report_path,reports);outputs.extend([str(p),str(report_path)])
            else:
                for arr in buffers.values():arr.flush()
                del buffers;products.append(product);outputs.extend(product[k] for k in ('data','delta','rfi'))
            rows+=n;print(stage,key[0],key[1],stem,n,flush=True)
    result=dict(scan=key[0],beam=key[1],rows=rows,outputs=outputs,products=products,stage=stage,fingerprint=fingerprint(stage,ds,params))
    dump(done,result);return result


def run_batch(recipe_path):
    recipe=read(recipe_path);root=Path(recipe['root']);root.mkdir(parents=True,exist_ok=True)
    ds=PreviewDataset(recipe['config'],recipe['selection'],memory_mib=128);params=recipe['params'];stage=recipe['stage']
    if ds.catalog_sha256!=recipe['catalog_sha256'] or source_hash()!=recipe['source_sha256']:
        raise ValueError('Source/catalog changed since preview; reopen the app before running')
    out=stage_path(root,stage,ds,params);out.mkdir(parents=True,exist_ok=True)
    with (root/'.batch.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        dump(out/'recipe.json',recipe)
        catalog=load_catalog(ds.cfg.get('rfi_catalog','packaged'),ds.config_path.parent);dump(out/'rfi_catalog_snapshot.json',catalog)
        started=time.time()
        try:
            if stage=='cube':
                previous=stage_path(root,'resample',ds,params)
                if not (previous/'summary.json').exists():raise ValueError('Run all resampling first with these same parameters')
                selection=dict(ds.selection,output_freq_ghz=output_grid(dict(ds.cfg,**params)).tolist())
                dump(out/'selection.json',selection);dump(out/'processed.json',read(previous/'processed.json'))
                cfg=dict(ds.cfg,**params,config_sha256=fingerprint(stage,ds,params))
                if cfg.get('baseline_components') is not None:cfg['sensitivity_missing_diagnostic']=True
                dump(root/'status.json',dict(stage=stage,state='running',phase='scan gridding',output=str(out)))
                grid(cfg,out)
                dump(root/'status.json',dict(stage=stage,state='running',phase='combination',output=str(out)))
                result=combine(cfg,out)
                if recipe.get('validate',False):
                    dump(root/'status.json',dict(stage=stage,state='running',phase='validation',output=str(out)))
                    result['validation']=validate(cfg,out)
                dump(out/'summary.json',result)
            else:
                dependency={'calibration':'reference','baseline':'calibration','resample':'baseline'}.get(stage)
                if dependency and not (stage_path(root,dependency,ds,params)/'summary.json').exists():
                    raise ValueError('Run all '+dependency+' first with these same parameters')
                keys=list(_groups(ds));results=[]
                def completed(result,number,total):
                    dump(root/'status.json',dict(stage=stage,state='running',completed=number,total=total,output=str(out),scheduler=str(out/'scheduler.json')))
                results=memory_map(run_group,keys,ds.cfg,stage,
                    worker_estimate(ds.cfg,stage,channels=max(ctx['entry']['channels'] for ctx in ds.files)),
                    initializer=init_worker,initargs=(str(ds.config_path),str(ds.selection_path),params,str(root),stage),
                    status_path=out/'scheduler.json',on_result=completed)
                results.sort(key=lambda x:(x['scan'],x['beam']))
                result=dict(stage=stage,scan_beams=len(results),spectra=sum(x['rows'] for x in results),seconds=time.time()-started,fingerprint=fingerprint(stage,ds,params))
                if stage=='resample':dump(out/'processed.json',results)
                dump(out/'summary.json',result)
            dump(root/'status.json',dict(stage=stage,state='complete',output=str(out),summary=result))
            print('COMPLETE',json.dumps(result),flush=True)
        except Exception as e:
            dump(root/'status.json',dict(stage=stage,state='failed',error=str(e),output=str(out)))
            raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--recipe',required=True);args=p.parse_args();run_batch(args.recipe)
