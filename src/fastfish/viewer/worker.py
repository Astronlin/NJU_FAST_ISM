"""Persistent SSH stdio worker. No Qt dependency; raw FITS remain remote."""
import argparse,base64,contextlib,io,json,sys,traceback,zipfile
from pathlib import Path
import numpy as np
from fastfish.processing.session import read_session,write_session,load_session_catalog
from fastfish.processing.preview import PreviewDataset,C,output_grid,resample
from fastfish.commands.shell import validate_params
from fastfish.rfi.catalog import load_catalog


def preview_arrays(ds,s,view=None):
    """Compute only the visible stage; an omitted view retains the old protocol."""
    if view not in (None,'Tsys calibration','Apply calibration','Baseline','Table'):
        raise ValueError('Unknown preview view: '+str(view))
    params=s['params'];index=s['index'];fi,ch,row=ds.address(index)
    ctx=ds.files[fi];local=ch*ds.cfg['chunk_rows']+row
    selected=dict(scan=ctx['scan'],beam=ctx['beam'],raw_row=int(ctx['coords']['ids'][local]),file=Path(ctx['entry']['path']).name)
    if view=='Tsys calibration':
        return ds.reference_display(fi,ch,params),dict(preview_view=view,selection_info=selected)
    if view=='Apply calibration':
        cal=ds.calibrated(fi,ch,params);sc=cal['scales'][row];f=cal['frequency_ghz']
        intervals=np.asarray(cal['intervals'],dtype=float).reshape(-1,2)
        return dict(velocity=(1-f*sc/ds.cfg['rest_freq_ghz'])*C,
            raw_pol=cal['raw_pol'][row],pol=cal['pol'][row],
            rfi_velocity=(1-intervals/1000*sc/ds.cfg['rest_freq_ghz'])*C),dict(preview_view=view,selection_info=selected)
    point=ds.point(index,params,selected_only=True);cal=point['cal'];fit=point['fit'];f=cal['frequency_ghz']
    arrays=dict(velocity=(1-f*cal['scales'][0]/ds.cfg['rest_freq_ghz'])*C,
        pol=cal['pol'][0],broad=fit['broad_pol'][0],local=fit['local_pol'][0],residual=fit['residual_pol'][0])
    if 'component_models' in fit:
        arrays['component_models']=np.stack([m[0] for m in fit['component_models']]) if fit['component_models'] else np.empty((0,len(f),2))
    if view in (None,'Table'):
        cfg=dict(ds.cfg,**params);q=output_grid(cfg)
        values,_=resample(fit['residual'],f,cal['scales'],cal['intervals'],cfg,q)
        arrays.update(output_velocity=(1-q/ds.cfg['rest_freq_ghz'])*C,output=values[0])
    if view is None:
        arrays.update(raw_pol=cal['raw_pol'][0],**ds.reference_display(fi,ch,params))
    intervals=np.asarray(cal['intervals'],dtype=float).reshape(-1,2)
    arrays['rfi_velocity']=(1-intervals/1000*cal['scales'][0]/ds.cfg['rest_freq_ghz'])*C
    return arrays,dict(fit_report=fit['report'],preview_view=view,selection_info=selected)


def encode_arrays(arrays,compress=True):
    """Lossless NPZ: avoid compression for local pipes, use fast deflate for SSH."""
    buffer=io.BytesIO()
    if not compress:
        np.savez(buffer,**arrays)
    else:
        with zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as archive:
            for name,array in arrays.items():
                with archive.open(name+'.npy','w',force_zip64=True) as member:
                    np.lib.format.write_array(member,np.asanyarray(array),allow_pickle=False)
    return buffer.getvalue()


def serve(path):
    ds=None;key=None
    for line in sys.stdin:
        try:
            req=json.loads(line);s=read_session(path)
            with contextlib.redirect_stdout(sys.stderr):
                if req['op']=='catalog':
                    edit=req.get('updates',{});s=load_session_catalog(path,edit['catalog_path'],edit.get('catalog_profile'),req['revision'])
                k=(s['config'],s['selection'])
                if key!=k:ds=PreviewDataset(*k,memory_mib=128);key=k
                if req.get('updates') and req['op']!='catalog':
                    updates=req['updates'];patch={}
                    if 'params' in updates:patch['params']=validate_params(ds.cfg,dict(s['params'],**updates['params']))
                    if 'index' in updates:
                        index=int(updates['index'])
                        if not 0<=index<len(ds.ra):raise ValueError('Spectrum outside selection')
                        patch['index']=index
                    s=write_session(path,patch,req['revision'])
                s['params']=validate_params(ds.cfg,s['params'])
                arrays={}
                info=dict(s.get('rfi_catalog_info',{}))
                if not info:
                    catalog=load_catalog(ds.cfg['rfi_catalog'])
                    info=dict(source=ds.cfg['rfi_catalog'],version=catalog['catalog_version'],profile=ds.cfg.get('rfi_profile'))
                s=dict(s,rfi_catalog_info=info,session_path=str(Path(path).resolve()))
                if req['op']=='map' or (req['op']=='catalog' and req.get('view')=='Footprint'):
                    names=list(ds.cfg['scans']);lookup={name:i for i,name in enumerate(names)}
                    groups=np.concatenate([np.full(len(ctx['coords']['ids']),lookup[ctx['scan']],dtype=np.min_scalar_type(max(0,len(names)-1))) for ctx in ds.files])
                    arrays=dict(ra=ds.ra,dec=ds.dec,scan_group=groups);s=dict(s,scan_names=names)
                elif req['op']=='region':
                    ids=ds.polygon(req['vertices']) if 'vertices' in req else ds.rectangle(*req['bounds']);r=ds.region_spectrum(ids,s['params'],max_spectra=128);arrays=dict(region_velocity=r['velocity'],region_mean=r['mean']);s=dict(s,region_used=r['used'],region_total=r['total'])
                elif req['op'] in ('point','catalog'):
                    arrays,extra=preview_arrays(ds,s,req.get('view'));s=dict(s,**extra)
                elif req['op']!='state':raise ValueError('Unknown operation')
            blob=encode_arrays(arrays,compress=req.get('compress',True))
            result=dict(state=s,arrays=base64.b64encode(blob).decode())
        except Exception as e:result=dict(error=str(e));traceback.print_exc(file=sys.stderr)
        print(json.dumps(result),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--session-state',required=True);a=p.parse_args();serve(a.session_state)
if __name__=='__main__':main()
