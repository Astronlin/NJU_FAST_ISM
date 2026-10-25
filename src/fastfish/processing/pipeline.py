"""Standalone, restartable steps for the agreed FAST Pulsar+Spectral workflow.

Only native SCI, extracted Tsys references and telemetry are inputs. Earlier
test cubes or processed spectra are never required. Science remains at native
resolution until baseline subtraction and the single final resampling.
"""
from fastfish._paths import PACKAGE_ROOT, DATA_ROOT
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import warnings
from fastfish.processing.scheduler import memory_map,worker_estimate,safe_block

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.utils import iers
from astropy.wcs import WCS

from fastfish.io.fits_stream import fixed_table, spectral_metadata,contiguous_rows
from fastfish.pointing.coordinates import reconstruct_pointing
from fastfish.pointing.doppler import frequency_scales
from fastfish.calibration.reference import reference_slice, telemetry_cache, smoothing_rows
from fastfish.processing.spectral import (frequency_gain, affine_overlap_rebin,
                                  gaussian_to_fwhm, affine_linear)
from fastfish.rfi.detection import (C, learn_priors, track_intervals, frequency_mask,
                           fit_linear_anchors)
from fastfish.baseline.broad import fit_global_broad
from fastfish.baseline.windows import normalize_windows
from fastfish.processing.preview import fit_polarization_baselines, output_grid, local_polarization_mean, resample
from fastfish.calibration.extracted import read_reference, reference_bounds
from fastfish.rfi.catalog import load_catalog, catalog_path, profile_priors, validate_catalog

STEPS = ('inventory', 'pointing', 'rfi', 'reduce', 'grid', 'combine', 'validate')


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n')
    temp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def source_hash():
    h = hashlib.sha256()
    package = PACKAGE_ROOT
    for p in sorted(package.rglob('*')):
        if p.is_file() and not p.name.startswith('._') and p.suffix in ('.py', '.fits', '.json'):
            h.update(str(p.relative_to(package)).encode()); h.update(p.read_bytes())
    return h.hexdigest()


def configure(path):
    path = Path(path).resolve()
    cfg = read(path)
    cfg['baseline_exclusions_kms']=normalize_windows(cfg.get('baseline_exclusions_kms',[]))
    cfg['n_t_apply']=smoothing_rows(cfg.get('n_t_apply',cfg.get('n_t',10)))
    cfg.pop('n_t',None)
    cfg['n_t_cal']=smoothing_rows(cfg.get('n_t_cal',10))
    if 'rfi_priors' in cfg:
        raise ValueError('Move inline rfi_priors into a shared catalog profile; use rfi_catalog and rfi_profile')
    reference_path = catalog_path(cfg.get('rfi_catalog', 'packaged'), path.parent)
    reference_bytes = reference_path.read_bytes()
    # Parse and hash the same bytes, then freeze them with this run.
    catalog = validate_catalog(json.loads(reference_bytes))
    cfg['rfi_catalog_path'] = str(reference_path)
    cfg['rfi_catalog_sha256'] = hashlib.sha256(reference_bytes).hexdigest()
    cfg['rfi_catalog_version'] = catalog['catalog_version']
    if cfg.get('rfi_profile') is not None:
        profile_priors(catalog, cfg['rfi_profile'], cfg['scans'])
    cfg['output_dir'] = str((path.parent / cfg.get('output_dir', 'output')).resolve())
    cfg['source_sha256'] = source_hash()
    cfg['config_sha256'] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    if not np.isfinite([cfg['protect_kms'],cfg['local_fit_kms'],cfg['alternative_fit_kms']]).all() or not 0<=cfg['protect_kms']<min(cfg['local_fit_kms'],cfg['alternative_fit_kms']):
        raise ValueError('Baseline fitting windows must exceed the nonnegative protected half-width')
    if not 1 <= cfg['chunk_rows'] <= 32 or cfg['workers'] < 1:
        raise ValueError('Use 1–32 native rows per chunk and at least one worker')
    if cfg['spectral_mode'] not in ('rebin', 'gaussian') or cfg['n_b'] < 1:
        raise ValueError('Invalid spectral mode or bin width')
    if cfg['native_fwhm_chan'] <= 0 or cfg['n_g'] < cfg['native_fwhm_chan']:
        raise ValueError('Gaussian target cannot be narrower than the assumed native LSF')
    root = Path(cfg['output_dir']); root.mkdir(parents=True, exist_ok=True)
    pinned = root / 'run_config.json'
    if pinned.exists() and read(pinned)['config_sha256'] != cfg['config_sha256']:
        raise ValueError('Configuration/source changed: choose a new output_dir')
    dump(pinned, cfg)
    (root / 'rfi_catalog_snapshot.json').write_bytes(reference_bytes)
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    iers.earth_orientation_table.set(iers.IERS_B(Table.read(DATA_ROOT / 'iers_b.fits')))
    return cfg, root


def inventory(cfg, root):
    items = cfg.get('inputs')
    if items is None:
        items = []
        for scan in cfg['scans']:
            for beam in cfg['beams']:
                folder = Path(cfg['input_root']) / scan / beam
                refs = list(folder.glob('Tsys_psr*.fits'))
                if len(refs) != 1:
                    raise ValueError('Expected one extracted Tsys reference: ' + str(folder))
                names = (folder / ('myspeclist_N_' + scan + '_' + beam)).read_text().splitlines()
                items.append(dict(scan=scan, beam=beam, reference=str(refs[0]),
                                  files=[str((folder / x).resolve()) if not Path(x).is_absolute() else x for x in names if x.strip()]))
    records = []; bounds = {s: [] for s in cfg['scans']}
    for item in items:
        if item['scan'] not in cfg['scans']:
            raise ValueError('Unexpected scan in inventory')
        if not Path(item['reference']).is_file():
            if not cfg.get('extracted_cal_root') or not Path(cfg.get('tcal_path','')).is_file():
                raise FileNotFoundError(item['reference'])
        entries = [spectral_metadata(p,time_reference=cfg.get('science_timestamp_reference','start')) for p in item['files']]
        if not entries:
            raise ValueError('Empty SCI list')
        for entry in entries:
            if entry['channels'] % 128 or not np.isclose(entry['native_df_hz'], cfg['native_df_hz'], rtol=1e-8):
                raise ValueError('Unexpected native frequency grid')
            bounds[item['scan']].extend(entry['mjd_first_last'])
        records.append(dict(scan=item['scan'], beam=item['beam'], reference=item['reference'], files=entries))
        print('INVENTORY', item['scan'], item['beam'], len(entries), flush=True)
    dump(root / 'inventory.json', dict(items=records, bounds=bounds, raw_psr_required=False))
    return dict(scan_beams=len(records), sci_files=sum(len(x['files']) for x in records))


def map_coordinates(ra, dec, frame='icrs'):
    """Transform only map positions; calibration/Doppler retain ICRS inputs."""
    if frame == 'icrs':return np.asarray(ra),np.asarray(dec)
    if frame != 'galactic':raise ValueError('coordinate_frame must be icrs or galactic')
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    position=SkyCoord(ra=np.asarray(ra)*u.deg,dec=np.asarray(dec)*u.deg,frame='icrs').galactic
    return position.l.deg,position.b.deg


def _auto_header(ra, dec, pixel_arcmin, frame='icrs', margin_pixels=None):
    lon,lat=map_coordinates(ra,dec,frame)
    good=np.isfinite(lon)&np.isfinite(lat);lon=lon[good];lat=lat[good]
    if not len(lon) or not np.isfinite(pixel_arcmin) or pixel_arcmin<=0:raise ValueError('Invalid automatic cube coordinates/pixel size')
    # Minimum circular longitude interval, including fields across 0/360 degrees.
    ordered=np.sort(lon%360);gaps=np.diff(np.r_[ordered,ordered[0]+360]);start=ordered[(np.argmax(gaps)+1)%len(ordered)]
    unwrapped=start+(lon-start)%360;centre=float((unwrapped.min()+unwrapped.max())/2%360)
    hdr=dict(NAXIS=2,NAXIS1=2,NAXIS2=2,CTYPE1='GLON-TAN' if frame=='galactic' else 'RA---TAN',
        CTYPE2='GLAT-TAN' if frame=='galactic' else 'DEC--TAN',CUNIT1='deg',CUNIT2='deg',
        CDELT1=-pixel_arcmin/60,CDELT2=pixel_arcmin/60,CRPIX1=1.,CRPIX2=1.,CRVAL1=centre,CRVAL2=float((lat.min()+lat.max())/2))
    if frame=='icrs':hdr['RADESYS']='ICRS'
    x,y=WCS(hdr).all_world2pix(lon,lat,0)
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Footprint exceeds a single TAN projection')
    # Same Gaussian support as _gridder: sigma = 1.5 arcmin / 2.35482.
    margin=(3*1.5/np.sqrt(8*np.log(2))/pixel_arcmin) if margin_pixels is None else float(margin_pixels)
    if not np.isfinite(margin) or margin<0:raise ValueError('Invalid cube kernel margin')
    xmin=int(np.floor(x.min()-margin+.5));xmax=int(np.ceil(x.max()+margin-.5))
    ymin=int(np.floor(y.min()-margin+.5));ymax=int(np.ceil(y.max()+margin-.5))
    hdr.update(NAXIS1=xmax-xmin+1,NAXIS2=ymax-ymin+1,CRPIX1=1.-xmin,CRPIX2=1.-ymin)
    return hdr


def pointing(cfg, root):
    inv = read(root / 'inventory.json'); telemetry = {}
    for scan in cfg['scans']:
        limits = inv['bounds'][scan]
        info = telemetry_cache(scan, cfg['pointing'][scan], min(limits), max(limits), root / 'pointing_cache')
        with np.load(info['path']) as f:
            telemetry[scan] = pd.DataFrame({k: f[k] for k in f.files})
    selections = []; all_ra = []; all_dec = []; excluded = 0
    for item in inv['items']:
        contexts = []
        tlo, thi = reference_bounds(item,cfg)
        for entry in item['files']:
            _, t = spectral_metadata(entry['path'], read_times=True,time_reference=cfg.get('science_timestamp_reference','start'))
            ra, dec = reconstruct_pointing(t, telemetry[item['scan']], beam=int(item['beam'][1:]))
            good = np.isfinite(ra) & np.isfinite(dec) & (t >= tlo) & (t <= thi)
            if cfg.get('header'):
                x, y = WCS(cfg['header']).all_world2pix(*map_coordinates(ra,dec,cfg.get('coordinate_frame','icrs')), 0); m = cfg['spatial_margin_pixels']
                # Footprint boundaries are pixel edges; expand by the kernel support.
                good &= (x >= -.5-m) & (x < cfg['header']['NAXIS1']-.5+m) & (y >= -.5-m) & (y < cfg['header']['NAXIS2']-.5+m)
            ids = np.flatnonzero(good)
            # Explicit smoke-test option; never enabled in the distributed examples.
            if cfg.get('smoke_rows_per_file'):
                ids = ids[np.unique(np.linspace(0,len(ids)-1,min(len(ids),cfg['smoke_rows_per_file'])).astype(int))] if len(ids) else ids
            excluded += int(len(t)-len(ids))
            if not len(ids):
                continue
            folder = root / 'selection' / item['scan'] / item['beam']; folder.mkdir(parents=True, exist_ok=True)
            out = folder / (Path(entry['path']).stem + '.npz')
            np.savez_compressed(out, ids=ids, mjd=t[ids], ra=ra[ids], dec=dec[ids],
                                frequency_scale=frequency_scales(t[ids], ra[ids], dec[ids], mode='spectrum'))
            all_ra.append(ra[ids]); all_dec.append(dec[ids])
            contexts.append(dict(entry=entry, coordinates=str(out)))
        if contexts:
            selections.append(dict(scan=item['scan'], beam=item['beam'], reference=item['reference'], files=contexts))
    if not selections:
        raise ValueError('No supported spectra in the requested footprint/reference times')
    hdr = cfg.get('header') or _auto_header(np.concatenate(all_ra), np.concatenate(all_dec), cfg['pixel_arcmin'],cfg.get('coordinate_frame','icrs'))
    q = output_grid(cfg)
    dump(root / 'selection.json', dict(items=selections, header=hdr, output_freq_ghz=q.tolist(),
         spectra=sum(len(x) for x in all_ra), excluded_rows=excluded))
    return dict(scan_beams=len(selections), spectra=sum(len(x) for x in all_ra), channels=len(q), header=hdr)


def _frequency(entry):
    return entry['freq_first_ghz'] + np.arange(entry['channels'])*entry['native_df_hz']/1e9


def identify_rfi(cfg, root):
    selection = read(root / 'selection.json')
    snapshot = root / 'rfi_catalog_snapshot.json'
    if hashlib.sha256(snapshot.read_bytes()).hexdigest() != cfg['rfi_catalog_sha256']:
        raise ValueError('Run RFI catalog snapshot fingerprint mismatch')
    catalog = load_catalog(snapshot)
    if cfg.get('rfi_profile') is not None:
        priors = profile_priors(catalog, cfg['rfi_profile'], cfg['scans'])
        origin = 'Shared catalog profile: ' + cfg['rfi_profile']
    else:
        profiles = []; scales = []; labels = []; common = None
        for item in selection['items']:
            ctx = item['files'][len(item['files'])//2]; f = _frequency(ctx['entry'])
            if common is not None and not np.array_equal(common, f):
                raise ValueError('RFI prior learning requires a shared native grid')
            common = f
            with np.load(ctx['coordinates']) as coords:
                jj = np.unique(np.linspace(0, len(coords['ids'])-1, 8).astype(int))
                ids = coords['ids'][jj]; mjd = coords['mjd'][jj]; fac = coords['frequency_scale'][jj]
            pack=read_reference(item,mjd.min(),mjd.max(),f,cfg)
            from scipy.interpolate import interp1d
            ref=interp1d(pack['mjd'],pack['smoothed'],axis=0,bounds_error=True);cf=pack['frequency_ghz']
            with fixed_table(ctx['entry']['path']) as (table, _):
                raw = table['DATA'][ids, :, :2].astype(float)
            protected = [((1-700/C)*cfg['rest_freq_ghz'], (1+700/C)*cfg['rest_freq_ghz'])]
            gain, _ = frequency_gain(raw, f, ref(mjd), cf, excluded_ranges=protected, degree=cfg['gain_degree'])
            profiles.extend((raw/gain).reshape(len(ids), -1, 128, 2).mean(axis=2))
            scales.extend(fac); labels.extend([item['scan']+'/'+item['beam']]*len(ids))
        priors = learn_priors(np.asarray(profiles), common.reshape(-1, 128).mean(axis=1),
                              np.asarray(scales), labels, cfg['rest_freq_ghz'],
                              pattern=catalog['patterns']['repeated_broad_features'])
        origin = 'New candidate priors from eight calibrated samples per selected scan/beam; review rfi_priors.json before reduction'
    for scan in cfg['scans']:
        if scan not in priors:
            raise ValueError('Missing RFI priors for ' + scan)
    dump(root / 'rfi_priors.json', priors)
    return dict(origin=origin, catalog_version=cfg['rfi_catalog_version'],
                catalog_sha256=cfg['rfi_catalog_sha256'],
                candidates_per_scan={s:len(p) for s,p in priors.items()})


def _reduce_item(job):
    cfg, root, item, q, priors = job; root = Path(root)
    folder = root / 'processed' / item['scan'] / item['beam']; folder.mkdir(parents=True, exist_ok=True)
    done = folder / 'done.json'
    if done.exists():
        record = read(done)
        if record['config_sha256'] != cfg['config_sha256'] or record.get('rfi_priors_sha256') != cfg['rfi_priors_sha256']:
            raise ValueError('Stale processed product')
        if all(Path(p).is_file() for z in record['products'] for p in z.values()):
            return record
    f = _frequency(item['files'][0]['entry']); ff = f.reshape(-1,128).mean(axis=1)
    times = []
    for ctx in item['files']:
        with np.load(ctx['coordinates']) as c:
            times.extend([float(c['mjd'][0]), float(c['mjd'][-1])])
    pack=read_reference(item,min(times),max(times),f,cfg)
    from scipy.interpolate import interp1d
    ref=interp1d(pack['mjd'],pack['smoothed'],axis=0,bounds_error=True);cf=pack['frequency_ghz']
    products = []; marked = 0; rows = 0; interval_records = []
    for ctx in item['files']:
        if not np.array_equal(f, _frequency(ctx['entry'])):
            raise ValueError('Native grid changed within beam')
        with np.load(ctx['coordinates']) as c:
            coords = {k:c[k] for k in c.files}
        n = len(coords['ids']); stem = Path(ctx['entry']['path']).stem
        files = {k:str(folder/(stem+'_'+k+'.npy')) for k in ('data','rfi','delta')}
        files['coordinates'] = ctx['coordinates']
        buffers = {k:np.lib.format.open_memmap(files[k],mode='w+',dtype='float32',shape=(n,len(q))) for k in ('data','rfi','delta')}
        streaming=contiguous_rows(coords['ids'])
        with fixed_table(ctx['entry']['path'],random_access=not streaming,sequential_access=streaming) as (tab, _):
            for k in range(0,n,cfg['chunk_rows']):
                sl = slice(k,min(k+cfg['chunk_rows'],n)); fac = coords['frequency_scale'][sl]
                raw = tab['DATA'][coords['ids'][sl],:,:2].astype(float)
                anchor = raw.mean(axis=(0,2)).reshape(-1,128).mean(axis=1)
                intervals, tracks = track_intervals(anchor,ff,priors[item['scan']],float(np.median(fac)),cfg['rest_freq_ghz'])
                masks = [((1-700/C)*cfg['rest_freq_ghz'],(1+700/C)*cfg['rest_freq_ghz'])]+[(a/1000,b/1000) for a,b in intervals]
                gain, _ = frequency_gain(raw,f,ref(coords['mjd'][sl]),cf,excluded_ranges=masks,degree=cfg['gain_degree'])
                calibrated_pol = raw/gain
                calibrated = calibrated_pol.mean(axis=2)
                vv = (1-f[None,:]*fac[:,None]/cfg['rest_freq_ghz'])*C
                if cfg['baseline']:
                    fitted = fit_polarization_baselines(calibrated_pol,f,fac,intervals,cfg,sensitivity=True)
                    residual=fitted['residual'];report=fitted['report']
                    coeff=fitted['local_coeff'];alt=fitted['alternative_coeff']
                else:
                    residual = calibrated; report = dict(baseline=False)
                    coeff = alt = np.zeros((len(fac),2,2))
                flag = frequency_mask(f,intervals,half_channel=cfg['native_df_hz']/2e6).astype(float)
                quality = np.broadcast_to(flag,residual.shape)
                if cfg['spectral_mode']=='rebin':
                    values = affine_overlap_rebin(residual,f,q,fac,axis=1,min_coverage=1.)[0]
                    fraction = affine_overlap_rebin(quality,f,q,fac,axis=1,min_coverage=1.)[0]
                else:
                    values = affine_linear(gaussian_to_fwhm(residual,cfg['n_g'],cfg['native_fwhm_chan']),f,q,fac,axis=1)
                    # Propagate marks through the same Gaussian and interpolation response.
                    fraction = affine_linear(gaussian_to_fwhm(quality,cfg['n_g'],cfg['native_fwhm_chan']),f,q,fac,axis=1)
                if not np.isfinite(values).all() or not np.isfinite(fraction).all():
                    raise ValueError('Unsupported/nonfinite output: '+stem)
                vout=(1-q/cfg['rest_freq_ghz'])*C
                buffers['data'][sl]=values; buffers['rfi'][sl]=np.clip(fraction,0,1)
                delta_native=local_polarization_mean(coeff,vv,cfg['local_fit_kms'])-local_polarization_mean(alt,vv,cfg['alternative_fit_kms'])
                buffers['delta'][sl]=resample(delta_native,f,fac,intervals,cfg,q)[0]
                marked += int((fraction>0).sum())
                interval_records.append(dict(input=ctx['entry']['path'],rows=[k,sl.stop],intervals_MHz=intervals,tracks=tracks,global_fit=report))
        for array in buffers.values():
            array.flush()
        del buffers
        rows += n; products.append(files)
        print('REDUCED',item['scan'],item['beam'],stem,n,flush=True)
    import gzip
    with gzip.open(str(folder/'rfi_intervals.json.gz'),'wt') as out:
        json.dump(interval_records,out)
    record=dict(scan=item['scan'],beam=item['beam'],products=products,rows=rows,marked_output_samples=marked,
                config_sha256=cfg['config_sha256'],rfi_priors_sha256=cfg['rfi_priors_sha256'])
    dump(done,record)
    return record


def reduce(cfg, root):
    selection=read(root/'selection.json'); priors=read(root/'rfi_priors.json'); q=np.asarray(selection['output_freq_ghz'])
    cfg=dict(cfg,rfi_priors_sha256=hashlib.sha256((root/'rfi_priors.json').read_bytes()).hexdigest())
    jobs=((cfg,str(root),item,q,priors) for item in selection['items'])
    records=memory_map(_reduce_item,jobs,cfg,'reduce',worker_estimate(cfg,'reduce'),status_path=root/'scheduler_reduce.json')
    dump(root/'processed.json',records)
    return dict(scan_beams=len(records),spectra=sum(r['rows'] for r in records))


def _gridder(header):
    import cygrid
    grid=cygrid.WcsGrid(header); sigma=(3/60)/2/np.sqrt(8*np.log(2))
    grid.set_kernel('gauss1d',(sigma,),3*sigma,sigma/2)
    grid.set_num_threads(1)
    return grid


def _grid_fields(cfg):
    return ('data','delta','rfi','sensitivity_missing') if cfg.get('sensitivity_missing_diagnostic',False) else ('data','delta','rfi')


def _grid_block(task):
    cfg,root,scan,a,b=task;root=Path(root)
    selection=read(root/'selection.json');records=read(root/'processed.json');h=selection['header']
    folder=root/'scan_grids'/scan;width=b-a;g=_gridder(h)
    for record in records:
        if record['scan']!=scan:continue
        for product in record['products']:
            arrays={name:np.load(product[name],mmap_mode='r') for name in _grid_fields(cfg)}
            with np.load(product['coordinates']) as c:
                for i in range(0,len(c['ra']),cfg['grid_row_block']):
                    sl=slice(i,i+cfg['grid_row_block'])
                    values=np.ascontiguousarray(np.concatenate([arrays[name][sl,a:b] for name in _grid_fields(cfg)],axis=1),dtype='float32')
                    g.grid(*map_coordinates(c['ra'][sl],c['dec'][sl],cfg.get('coordinate_frame','icrs')),values)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning);data=g.get_datacube()
    weights=g.get_weights()
    if not np.all(weights==weights[0]):raise ValueError('Retained-channel geometric weights unexpectedly differ')
    digest=hashlib.sha256(np.ascontiguousarray(weights[0]).tobytes()).hexdigest()
    if a==0:np.save(str(folder/'weights.npy'),weights[0])
    for j,name in enumerate(_grid_fields(cfg)):
        output=np.load(str(folder/(name+'.npy')),mmap_mode='r+')
        output[a:b]=data[j*width:(j+1)*width];output.flush();del output
    return dict(scan=scan,start=a,stop=b,weights_sha256=digest)


def grid(cfg, root):
    root=Path(root);selection=read(root/'selection.json');records=read(root/'processed.json');h=selection['header']
    nq=len(selection['output_freq_ghz']);shape=(nq,h['NAXIS2'],h['NAXIS1']);partials=[];tasks=[]
    width=safe_block(cfg,cfg['grid_channel_block'],h['NAXIS1']*h['NAXIS2'],factor=160 if cfg.get('sensitivity_missing_diagnostic',False) else 120)
    for scan in sorted(set(r['scan'] for r in records)):
        folder=root/'scan_grids'/scan;folder.mkdir(parents=True,exist_ok=True)
        paths={name:str(folder/(name+'.npy')) for name in (*_grid_fields(cfg),'weights')}
        for name in _grid_fields(cfg):
            output=np.lib.format.open_memmap(paths[name],mode='w+',dtype='float32',shape=shape);del output
        partials.append(dict(scan=scan,**paths))
        tasks.extend((cfg,str(root),scan,a,min(nq,a+width)) for a in range(0,nq,width))
    def completed(result,number,total):
        dump(root/'grid_progress.json',dict(completed_blocks=number,total_blocks=total,scan=result['scan'],last_channels=[result['start'],result['stop']],channel_block=width))
    results=memory_map(_grid_block,tasks,cfg,'grid',int(worker_estimate(cfg,'grid',pixels=h['NAXIS1']*h['NAXIS2'],width=width)*(4/3 if cfg.get('sensitivity_missing_diagnostic',False) else 1)),status_path=root/'scheduler_grid.json',on_result=completed)
    for part in partials:
        hashes={r['weights_sha256'] for r in results if r['scan']==part['scan']}
        if len(hashes)!=1:raise ValueError('Spectral blocks have different weights')
        print('GRIDDED',part['scan'],flush=True)
    dump(root/'scan_grids.json',partials)
    return dict(scans=len(partials),shape=list(shape),channel_block=width)


def _header(cfg, selection):
    h=fits.Header(selection['header']); q=np.asarray(selection['output_freq_ghz'])
    for key,value in dict(CTYPE3='FREQ',CUNIT3='Hz',CRPIX3=1.,CRVAL3=q[0]*1e9,CDELT3=(q[1]-q[0])*1e9,
         RESTFRQ=cfg['rest_freq_ghz']*1e9,SPECSYS='LSRK',BUNIT='K',BASELINE=cfg['baseline'],
         BASEPOL='AA,BB',BASEMOD='SINE+LIN' if cfg['baseline'] else 'NONE',BASEVMIN=-float(cfg['local_fit_kms']),BASEVMAX=float(cfg['local_fit_kms']),PROTVMIN=-float(cfg['protect_kms']),PROTVMAX=float(cfg['protect_kms']),
         DOPMODE='SPECTRUM',RFIEXCIS=False,NSCANS=len(cfg['scans']),SPCMODE=cfg['spectral_mode'],
         NATDFHZ=cfg['native_df_hz'],REBINN=cfg['n_b'],GFWHMCH=cfg['n_g'],CFGHASH=cfg['config_sha256']).items():
        h[key]=value
    if cfg.get('standing_wave_correction',False):
        h['BASERIPL']=True
        h['RIPPMIN']=1.04;h['RIPPMAX']=1.14
        h['MINREFOK']=float(cfg.get('minimum_positive_reference_fraction',.9))
        h['QAREJECT']=int(selection.get('calibration_quality_excluded_rows',0))
        h.add_history('Conditional short-period ripple: >=10 percent held-out RMS improvement and >=0.03 K amplitude; protected H I excluded.')
        h.add_history('Rows failing calibration-reference positivity quality criterion excluded from this comparison only; original data retained.')
    if cfg.get('baseline_components') is not None:
        h['BASEMOD']='COMPONENTS'
        h['NBASECMP']=len(cfg['baseline_components']) if cfg.get('baseline',True) else 0
        for key in ('BASEVMIN','BASEVMAX'):
            if key in h:del h[key]
        h.add_history('Ordered configurable components; fit/subtract only within each component range. See saved recipe.json for exact functions and ranges.')
    h['SCITIME']='MIDPOINT'
    h['RAWTIME']=cfg.get('science_timestamp_reference','start').upper()
    h['TIMEOFF']='EXPOSURE/2' if cfg.get('science_timestamp_reference','start')=='start' else 'NONE'
    h.add_history('Science calibration, pointing and Doppler use integration midpoint times; raw timestamps unchanged.')
    h.add_history('RFI excluded from gain/baseline fits only; marked science values retained. Geometric kernel weights.')
    return h


def combine(cfg, root):
    selection=read(root/'selection.json'); partials=read(root/'scan_grids.json'); h=selection['header']
    shape=(len(selection['output_freq_ghz']),h['NAXIS2'],h['NAXIS1'])
    block_width=safe_block(cfg,cfg['grid_channel_block'],h['NAXIS1']*h['NAXIS2'],factor=48)
    total=np.zeros(shape[1:],dtype=float)
    for part in partials:total+=np.load(part['weights'],mmap_mode='r')
    fits.PrimaryHDU(total.astype('float32'),fits.Header(h)).writeto(root/'coverage.fits',overwrite=True)
    for field,name in [('data',cfg['name']+'_cube.fits'),('rfi','rfi_contribution_fraction.fits'),('delta','baseline_sensitivity.fits')]+([('sensitivity_missing','baseline_sensitivity_missing_fraction.fits')] if cfg.get('sensitivity_missing_diagnostic',False) else []):
        scratch=root/(field+'_combined.npy'); out=np.lib.format.open_memmap(scratch,mode='w+',dtype='float32',shape=shape)
        for a in range(0,shape[0],block_width):
            b=min(shape[0],a+block_width); num=np.zeros((b-a,)+shape[1:])
            for part in partials:
                w=np.load(part['weights'],mmap_mode='r'); d=np.load(part[field],mmap_mode='r')[a:b]
                num+=np.where(w>0,d,0)*w
            combined=np.divide(num,total,out=np.full_like(num,np.nan),where=total>0)
            if field=='delta' and cfg.get('sensitivity_missing_diagnostic',False):
                # Placeholder zeros NEVER become reported zero sensitivity.
                missing=np.zeros_like(num)
                for part in partials:
                    w=np.load(part['weights'],mmap_mode='r');d=np.load(part['sensitivity_missing'],mmap_mode='r')[a:b]
                    missing+=np.where(w>0,d,0)*w
                combined[missing>0]=np.nan
            out[a:b]=combined
        out.flush(); hh=_header(cfg,selection)
        if field=='rfi':
            hh['BUNIT']='fraction'; hh.add_history('Weighted fraction marked RFI/candidate; not contamination probability or RFI amplitude.')
        elif field=='sensitivity_missing':
            hh['BUNIT']='fraction';hh.add_history('Geometric weighted fraction with unavailable alternative baseline fit; sensitivity is NaN wherever this is positive.')
        elif field=='delta':
            if cfg.get('baseline_components') is None:hh.add_history('Local +/-'+str(cfg['alternative_fit_kms'])+' minus +/-'+str(cfg['local_fit_kms'])+' fit; same global model/masks. Not a statistical uncertainty.')
        if field=='delta' and cfg.get('sensitivity_missing_diagnostic',False):hh.add_history('Alternative-window sensitivity is undefined for configurable components; see missing-fraction map.' if cfg.get('baseline_components') is not None else 'NaN where any contributing row lacks eight clean alternative-fit anchors on either side; see baseline_sensitivity_missing_fraction.fits.')
        fits.PrimaryHDU(out,hh).writeto(root/name,overwrite=True)
        del out; scratch.unlink()
    result=dict(shape=list(shape),spectra=selection['spectra'],scans=len(partials),config_sha256=cfg['config_sha256'])
    dump(root/'summary.json',result)
    return result


def validate(cfg, root):
    selection=read(root/'selection.json'); h=selection['header']; q=np.asarray(selection['output_freq_ghz'])
    indices=np.unique(np.linspace(0,len(q)-1,7).astype(int)); g=_gridder(h)
    for record in read(root/'processed.json'):
        for p in record['products']:
            data=np.load(p['data'],mmap_mode='r')
            with np.load(p['coordinates']) as c:
                for i in range(0,len(c['ra']),cfg['grid_row_block']):
                    sl=slice(i,i+cfg['grid_row_block'])
                    g.grid(*map_coordinates(c['ra'][sl],c['dec'][sl],cfg.get('coordinate_frame','icrs')),np.ascontiguousarray(data[sl][:,indices]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning); direct=g.get_datacube()
    with fits.open(root/(cfg['name']+'_cube.fits'),memmap=True) as cube:
        merged=cube[0].data[indices]; valid=np.isfinite(direct)&np.isfinite(merged)
        if not np.array_equal(np.isfinite(direct),np.isfinite(merged)):
            raise ValueError('Direct/merged finite-support mismatch')
        diff=direct[valid]-merged[valid]
        if not len(diff):
            raise ValueError('No covered map pixels in the selected sample')
        if np.max(abs(diff))>1e-3:
            raise ValueError('Direct/merged arithmetic check failed')
        finite=0
        for a in range(0,len(q),cfg['grid_channel_block']):
            block=cube[0].data[a:a+cfg['grid_channel_block']]
            if np.isinf(block).any():
                raise ValueError('Infinite cube values')
            finite+=int(np.isfinite(block).sum())
        fraction=finite/cube[0].data.size
    with fits.open(root/'rfi_contribution_fraction.fits',memmap=True) as marked:
        for a in range(0,len(q),cfg['grid_channel_block']):
            block=marked[0].data[a:a+cfg['grid_channel_block']]; x=block[np.isfinite(block)]
            if len(x) and (x.min() < -1e-6 or x.max() > 1+1e-6):
                raise ValueError('Invalid RFI fraction')
    sensitivity_unavailable=0
    if cfg.get('sensitivity_missing_diagnostic',False):
        with fits.open(root/'baseline_sensitivity_missing_fraction.fits',memmap=True) as marked, fits.open(root/'baseline_sensitivity.fits',memmap=True) as sensitivity:
            for a in range(0,len(q),cfg['grid_channel_block']):
                missing=marked[0].data[a:a+cfg['grid_channel_block']];delta=sensitivity[0].data[a:a+cfg['grid_channel_block']]
                finite=np.isfinite(missing);x=missing[finite]
                if len(x) and (x.min() < -1e-6 or x.max()>1+1e-6):raise ValueError('Invalid sensitivity-missing fraction')
                if not np.array_equal(np.isfinite(delta),finite&(missing==0)):raise ValueError('Sensitivity missingness/support mismatch')
                sensitivity_unavailable+=int((finite&(missing>0)).sum())
    report=dict(passed=True,direct_merge_max_K=float(np.max(abs(diff))),direct_merge_RMS_K=float(np.sqrt(np.mean(diff**2))),
                sampled_channels=indices.tolist(),finite_voxel_fraction=fraction,sensitivity_unavailable_voxels=sensitivity_unavailable,
                limitations='Arithmetic checks do not establish baseline truth, absolute calibration, or RFI-free H I wings.')
    dump(root/'validation.json',report)
    return report


def run_stage(config, stage):
    if stage not in STEPS:
        raise ValueError('Unknown step: '+stage)
    cfg, root=configure(config)
    with (root/'.pipeline.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        result=globals()['identify_rfi' if stage=='rfi' else stage](cfg,root)
        dump(root/(stage+'_report.json'),result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=STEPS+('all',))
    parser.add_argument('--config',required=True)
    args=parser.parse_args()
    for stage in STEPS if args.stage=='all' else (args.stage,):
        print(stage,json.dumps(run_stage(args.config,stage)),flush=True)


if __name__=='__main__':
    main()
