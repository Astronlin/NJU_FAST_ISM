"""Auditable downstream comparison using frozen calibration products.

Run --stage baseline, resample, then cube. The source run is immutable.
This comparison requires a saved reference-quality audit and records exclusions
separately from RFI. No calibration reference is regenerated or interpolated
across invalid intervals to invent replacement science values.
"""
from pathlib import Path
import argparse,fcntl,hashlib,json,time
import numpy as np
from fastfish.processing.pipeline import read,dump,source_hash,grid,combine,validate
from fastfish.processing.preview import (PreviewDataset,reconstruct_calibrated,fit_polarization_baselines,
    reconstruct_polarization_baselines,local_linear_correction,resample,output_grid)
from fastfish.baseline.standing_wave import fit_residual_wave,evaluate_wave
from fastfish.io.fits_stream import fixed_table,contiguous_rows
from fastfish.rfi.detection import track_intervals,fit_linear_anchors,C
from fastfish.processing.scheduler import memory_map,worker_estimate

_WORK=None

def _stage(source,name):
    matches=[p.parent for p in (source/'stages'/name).glob('*/summary.json')]
    if len(matches)!=1:raise ValueError('Expected one unambiguous completed source stage: '+name)
    return matches[0]


def initialize(source,root,stage):
    global _WORK
    source=Path(source);root=Path(root);ds=PreviewDataset(source/'config.json',source/'selection.json',128)
    groups={}
    for i,c in enumerate(ds.files):groups.setdefault((c['scan'],c['beam']),[]).append(i)
    _WORK=(ds,source,root,stage,groups,_stage(source,'calibration'))


def _intervals(raw,f,scales,priors,cfg):
    intervals=[];records=[]
    for p in range(2):
        anchor=raw[:,:,p].mean(axis=0).reshape(-1,128).mean(axis=1)
        mask,tracks=track_intervals(anchor,f.reshape(-1,128).mean(axis=1),priors,float(np.median(scales)),cfg['rest_freq_ghz'],protect=cfg['protect_kms'],catalog=cfg.get('rfi_catalog','packaged'))
        intervals.extend(mask);records.append(dict(polarization=['AA','BB'][p],tracks=tracks))
    merged=[]
    for lo,hi in sorted(intervals):
        if merged and lo<=merged[-1][1]:merged[-1][1]=max(hi,merged[-1][1])
        else:merged.append([lo,hi])
    return merged,records


def process_group(key):
    ds,source,root,stage,groups,calroot=_WORK;cfg=ds.cfg;scan,beam=key
    folder=root/stage/scan/beam;folder.mkdir(parents=True,exist_ok=True);done=folder/'done.json'
    if done.exists():
        record=read(done)
        if all(Path(p).is_file() for p in record['outputs']):return record
    outputs=[];products=[];rows=0;excluded=0;detected=0
    for fi in groups[tuple(key)]:
        ctx=ds.files[fi];c=ctx['coords'];f=ds.frequency(fi);n=len(c['ids']);stem=Path(ctx['entry']['path']).stem
        with np.load(calroot/scan/beam/(stem+'_calibration.npz')) as z:coeff=z['coefficients'];edges=z['reference_edges_ghz'];factors=z['frequency_scale']
        if len(coeff)!=n:raise ValueError('Frozen calibration/SCI row mismatch')
        basepath=root/'baseline'/scan/beam/(stem+'_baseline.npz');maskpath=basepath.with_name(stem+'_intervals.json')
        if stage=='baseline':
            quality=source/'baseline_review/reference_quality'/(stem+'_reference_quality.npz')
            with np.load(quality) as z:
                if not np.array_equal(z['mjd'],c['mjd']):raise ValueError('Reference quality times mismatch')
                fraction=z['nonpositive_reference_fraction'].copy()
            valid=np.isfinite(fraction)&(fraction<=.1)
            gc=np.full((n,2,3),np.nan);lc=np.full((n,2,2),np.nan);alt=np.full_like(lc,np.nan);wave=np.zeros_like(lc)
            periods=[];waveperiods=[];reports=[];masks=[]
        else:
            with np.load(basepath) as z:base={k:z[k] for k in z.files}
            valid=base['row_valid'];masks=read(maskpath);nkeep=int(valid.sum());q=output_grid(cfg)
            if not nkeep:
                excluded+=n;continue
            product={name:str(folder/(stem+'_'+name+'.npy')) for name in ['data','delta','rfi','sensitivity_missing']};product['coordinates']=str(folder/(stem+'_coordinates.npz'))
            coordinates={k:(v[valid] if np.ndim(v)>0 and len(v)==n else v) for k,v in c.items()};np.savez_compressed(product['coordinates'],**coordinates)
            buffers={name:np.lib.format.open_memmap(product[name],mode='w+',dtype='float32',shape=(nkeep,len(q))) for name in ['data','delta','rfi','sensitivity_missing']};cursor=0
        streaming=contiguous_rows(c['ids'])
        with fixed_table(ctx['entry']['path'],random_access=not streaming,sequential_access=streaming) as (tab,_):
            for chunk,a in enumerate(range(0,n,cfg['chunk_rows'])):
                b=min(n,a+cfg['chunk_rows']);ids=np.flatnonzero(valid[a:b])+a
                if not len(ids):
                    if stage=='baseline':periods.append([1.,1.]);waveperiods.append([1.09,1.09]);reports.append(dict(excluded_all_rows=True));masks.append(dict(intervals=[],tracks=[]))
                    continue
                raw=tab['DATA'][c['ids'][ids],:,:2].astype(float);fac=factors[ids];pol=reconstruct_calibrated(raw,f,coeff[ids],edges)
                if stage=='baseline':
                    intervals,tracks=_intervals(raw,f,fac,ds.priors[scan],cfg)
                    fit=fit_polarization_baselines(pol,f,fac,intervals,cfg,sensitivity=False)
                    gc[ids]=fit['global_coeff'];lc[ids]=fit['local_coeff'];wr=[];wp=[]
                    v=(1-f[None,:]*fac[:,None]/cfg['rest_freq_ghz'])*C
                    for p in range(2):
                        dl,wc,info=fit_residual_wave(fit['residual_pol'][:,:,p],f,fac,cfg['rest_freq_ghz'],intervals,limit=cfg['local_fit_kms'],protect=cfg['protect_kms'])
                        if info['detected']:
                            lc[ids,p]+=dl;wave[ids,p]=wc
                            detected+=1
                        model=evaluate_wave(wc,info['period_MHz'],f,v,cfg['local_fit_kms'])
                        alternative,alternative_info=fit_linear_anchors(pol[:,:,p]-fit['broad_pol'][:,:,p]-model,f,fac,cfg['rest_freq_ghz'],intervals,fit_limit=cfg['alternative_fit_kms'],protect=cfg['protect_kms'],extra_exclusions=cfg.get('baseline_exclusions_kms',[]),allow_missing=True)
                        alt[ids,p]=alternative
                        info['sensitivity_unavailable_rows']=int((~np.isfinite(alternative).all(axis=1)).sum())
                        info['alternative_minimum_left_anchors']=alternative_info['minimum_left_anchors']
                        info['alternative_minimum_right_anchors']=alternative_info['minimum_right_anchors']
                        wr.append(info);wp.append(info['period_MHz'])
                    periods.append([fit['report']['polarizations'][label]['period_MHz'] for label in ['AA','BB']]);waveperiods.append(wp);reports.append(dict(broad=fit['report'],ripple=wr));masks.append(dict(intervals=intervals,tracks=tracks))
                else:
                    residual,delta=reconstruct_polarization_baselines(pol,f,fac,base['global_coeff'][ids],base['local_coeff'][ids],base['period_MHz'][chunk],cfg,base['alternative_coeff'][ids])
                    v=(1-f[None,:]*fac[:,None]/cfg['rest_freq_ghz'])*C
                    for p in range(2):residual-=evaluate_wave(base['wave_coeff'][ids,p],base['wave_period_MHz'][chunk,p],f,v,cfg['local_fit_kms'])/2
                    intervals=masks[chunk]['intervals'];values,fraction=resample(residual,f,fac,intervals,cfg,q);sl=slice(cursor,cursor+len(ids))
                    buffers['data'][sl]=values;buffers['rfi'][sl]=fraction;missing=~np.isfinite(base['alternative_coeff'][ids]).all(axis=(1,2));delta[missing]=0.
                    if not np.isfinite(delta).all():raise ValueError('Nonfinite supported sensitivity')
                    buffers['delta'][sl]=resample(delta,f,fac,intervals,cfg,q)[0];buffers['sensitivity_missing'][sl]=missing[:,None];cursor+=len(ids)
        if stage=='baseline':
            np.savez_compressed(basepath,global_coeff=gc,local_coeff=lc,alternative_coeff=alt,period_MHz=periods,wave_coeff=wave,wave_period_MHz=waveperiods,row_valid=valid,nonpositive_reference_fraction=fraction,sensitivity_available=np.isfinite(alt).all(axis=(1,2)))
            dump(maskpath,masks);rp=folder/(stem+'_fit_reports.json');dump(rp,reports);outputs.extend(map(str,[basepath,maskpath,rp]));rows+=n
        else:
            assert cursor==nkeep
            for arr in buffers.values():arr.flush()
            del buffers;products.append(product);outputs.extend(product.values());rows+=nkeep
        excluded+=int((~valid).sum());print(stage,scan,beam,stem,'valid',int(valid.sum()),'excluded',int((~valid).sum()),flush=True)
    record=dict(scan=scan,beam=beam,rows=rows,excluded_rows=excluded,outputs=outputs,products=products,detected_ripple_polarization_chunks=detected);dump(done,record);return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-run',required=True);parser.add_argument('--run-root',required=True);parser.add_argument('--stage',choices=['baseline','resample','cube'],required=True);args=parser.parse_args()
    source=Path(args.source_run).resolve();root=Path(args.run_root).resolve()
    if source==root:raise ValueError('Correction requires a separate comparison workspace')
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.correction.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        cfg=read(source/'config.json');cfg.update(name=cfg.get('name',source.name)+'_corrected',workers=4,max_workers=4,io_workers=4,standing_wave_correction=True,minimum_positive_reference_fraction=.9,sensitivity_missing_diagnostic=True)
        provenance=dict(source_run=str(source),source_cube=str(_stage(source,'cube')),frozen_calibration=str(_stage(source,'calibration')),frozen_reference=str(_stage(source,'reference')),source_sha256=source_hash(),settings=cfg,policy='Conditional 1.04-1.14 MHz ripple; new observed-spectrum RFI detections; exclude rows with >10 percent nonpositive/nonfinite reference anchors from comparison gridding only. Original SCI/CAL/products preserved.')
        fingerprint=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest();cfg['config_sha256']=fingerprint
        if (root/'correction.json').exists() and read(root/'correction.json')!=provenance:raise ValueError('Correction settings/source changed; do not mix partial results')
        dump(root/'correction.json',provenance);dump(root/'config.json',cfg);stage=args.stage;folder=root/stage;folder.mkdir(exist_ok=True);started=time.time();dump(root/'status.json',dict(stage=stage,state='running'))
        try:
            if stage in ['baseline','resample']:
                if stage=='resample' and not (root/'baseline/summary.json').exists():raise ValueError('Complete baseline first')
                initialize(source,root,stage);ds=_WORK[0];keys=list(_WORK[4])
                def completed(record,number,total):dump(root/'status.json',dict(stage=stage,state='running',completed=number,total=total,output=str(folder)))
                result=memory_map(process_group,keys,cfg,stage,worker_estimate(cfg,stage,channels=max(c['entry']['channels'] for c in ds.files)),initializer=initialize,initargs=(str(source),str(root),stage),status_path=folder/'scheduler.json',on_result=completed)
                summary=dict(stage=stage,groups=len(result),rows=sum(x['rows'] for x in result),excluded_rows=sum(x['excluded_rows'] for x in result),ripple_polarization_chunks=sum(x['detected_ripple_polarization_chunks'] for x in result),seconds=time.time()-started)
                if stage=='resample':dump(folder/'processed.json',result)
            else:
                if not (root/'resample/summary.json').exists():raise ValueError('Complete resampling first')
                records=read(root/'resample/processed.json');sel=read(source/'selection.json');sel['spectra']=sum(r['rows'] for r in records);sel['calibration_quality_excluded_rows']=sum(r['excluded_rows'] for r in records)
                dump(folder/'selection.json',sel);dump(folder/'processed.json',records);grid(cfg,folder);dump(root/'status.json',dict(stage=stage,state='running',phase='combination'));summary=combine(cfg,folder);dump(root/'status.json',dict(stage=stage,state='running',phase='validation'));summary['validation']=validate(cfg,folder)
            dump(folder/'summary.json',summary);dump(root/'status.json',dict(stage=stage,state='complete',output=str(folder),summary=summary));print('COMPLETE',json.dumps(summary),flush=True)
        except Exception as error:
            dump(root/'status.json',dict(stage=stage,state='failed',error=str(error)));raise

if __name__=='__main__':main()
