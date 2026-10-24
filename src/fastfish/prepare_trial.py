"""Prepare a region's coordinate-only selection, optionally reusing audited telemetry results."""
from pathlib import Path
import argparse
import json
from collections import Counter
import numpy as np
from astropy.wcs import WCS
from .fits_stream import spectral_metadata
from .pipeline import dump,run_stage
from .interactive_core import output_grid


def prepare(config_path,root,seed=None):
    config_path=Path(config_path).resolve();root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    cfg=json.loads(config_path.read_text());allowed={x for item in cfg['inputs'] for x in item['files']}
    if seed is None:
        run_stage(config_path,'inventory');run_stage(config_path,'pointing')
        original=config_path.parent/cfg['output_dir']/'selection.json'
        dump(root/'selection.json',json.loads(original.read_text()));return
    old=json.loads(Path(seed).read_text());w=WCS(cfg['header']);m=cfg['spatial_margin_pixels'];items=[];counts=Counter({'0002':0,'0003':0,'0004':0});spectra=0
    nx=cfg['header']['NAXIS1'];ny=cfg['header']['NAXIS2']
    for item in old['selected']:
        contexts=[]
        for record in item['files']:
            path=record['entry']['path']
            if path not in allowed:continue
            number=Path(path).stem.split('_N_')[-1]
            if number not in ('0002','0003','0004'):raise ValueError('Unexpected file number: '+path)
            ids=np.asarray(record['selected_indices'])
            # Coordinate fields only: never import legacy calibrated science values.
            with np.load(record['previous_product']) as p:c={key:p[key][ids].copy() for key in ('ra','dec','mjd','frequency_scale')}
            x,y=w.all_world2pix(c['ra'],c['dec'],0);keep=(x>=-.5-m)&(x<nx-.5+m)&(y>=-.5-m)&(y<ny-.5+m)
            if not keep.any():continue
            c={k:v[keep] for k,v in c.items()};entry,t=spectral_metadata(path,read_times=True)
            rawids=np.searchsorted(t,c['mjd'])
            if not np.array_equal(t[rawids],c['mjd']):raise ValueError('Coordinate/native time mismatch')
            c['ids']=rawids;folder=root/'coordinates'/item['scan']/item['beam'];folder.mkdir(parents=True,exist_ok=True)
            out=folder/(Path(path).stem+'.npz');np.savez_compressed(out,**c)
            contexts.append(dict(entry=entry,coordinates=str(out)));spectra+=len(rawids);counts[number]+=len(rawids)
        if contexts:items.append(dict(scan=item['scan'],beam=item['beam'],reference=item['reference'],files=contexts))
    if not items:raise ValueError('Empty trial selection')
    selection=dict(items=items,header=cfg['header'],output_freq_ghz=output_grid(cfg).tolist(),spectra=spectra,file_number_counts=dict(counts))
    dump(root/'selection.json',selection)
    report=dict(spectra=spectra,scan_beams=len(items),SCI_files=sum(len(x['files']) for x in items),spectra_by_file_number=dict(counts),
        coordinate_seed=str(seed),scope='Reused audited RA/Dec/MJD/Doppler factors only. Science/calibration read from native inputs. Preview functions do not write.',
        size_arcmin=[81,78],kernel_margin_arcmin=m*cfg['pixel_arcmin'])
    dump(root/'selection_report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--root',required=True);p.add_argument('--coordinate-seed');a=p.parse_args();prepare(a.config,a.root,a.coordinate_seed)
