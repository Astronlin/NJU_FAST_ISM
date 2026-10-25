"""Bounded regeneration of missing Tsys references from existing extracted CAL.

No raw PSR extraction. This retains the legacy gaincal diode factor, configurable
ON-OFF filter, 16-channel Tcal filter, float32 stored-reference rounding, and
a configurable subsequent temperature-reference time filter (n_t_apply, default 10).
"""
from datetime import datetime,timezone
from functools import lru_cache
from pathlib import Path
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d
from astropy.io import fits
from fastfish.calibration.reference import reference_data,bracket,smoothing_rows
from fastfish.io.fits_stream import fixed_table


@lru_cache(maxsize=128)
def cal_index(root,scan,beam):
    result=[]
    for p in (Path(root)/scan).glob('*/*-'+beam+'_*.fits*'):
        match=re.search(r'_(\d{6})_(\d{6})(?:_cal)?\.fits(?:\.gz)?$',p.name)
        if match:
            stamp=datetime.strptime(''.join(match.groups()),'%y%m%d%H%M%S').replace(tzinfo=timezone.utc)
            result.append((stamp.timestamp()/86400+40587,str(p)))
    if not result:raise FileNotFoundError('No extracted CAL files for '+scan+'/'+beam)
    return tuple(sorted(result))


@lru_cache(maxsize=32)
def diode_reference(path,beam):
    with fits.open(path,memmap=True) as h:
        frequency=h[1].data['freq'][0].copy()/1e3
        values=h[1].data['tcal'][0,int(beam[1:])-1].astype(float).T
    return frequency,uniform_filter1d(values,16,axis=1,mode='nearest')


@lru_cache(maxsize=96)
def extracted_band(name,lower,upper):
    """Cache only the relevant ~260 PSR channels, never whole gzipped files."""
    with fits.open(name,memmap=False) as h:
        nf=int(h[0].header['OBSNCHAN']);df=float(h[1].header['CHAN_BW'])
        all_f=np.arange(float(h[0].header['OBSFREQ'])-float(h[0].header['OBSBW'])/2,
                        float(h[0].header['OBSFREQ'])+float(h[0].header['OBSBW'])/2,df)/1e3
        if len(all_f)!=nf:raise ValueError('Extracted CAL frequency metadata mismatch')
        a=max(0,np.searchsorted(all_f,lower)-2);b=min(nf,np.searchsorted(all_f,upper)+3)
        tab=h[2].data
        return tab['MJD'].copy(),tab['CAL'][:,1,:2,a:b].astype('float32'),tab['CAL'][:,0,:2,a:b].astype('float32'),all_f[a:b]


def rebuild_reference(paths,tmin,tmax,science_freq,tcal_path,beam,n_t=10,n_t_cal=10):
    n_t=smoothing_rows(n_t);n_t_cal=smoothing_rows(n_t_cal)
    times=[];ons=[];offs=[];cf=None
    for name in paths:
        tt,on,off,frequency=extracted_band(str(name),float(science_freq[0]),float(science_freq[-1]))
        if cf is not None and not np.array_equal(cf,frequency):raise ValueError('CAL frequency grids differ')
        cf=frequency;times.append(tt);ons.append(on);offs.append(off)
    tt=np.concatenate(times);on=np.concatenate(ons);off=np.concatenate(offs);order=np.argsort(tt)
    tt=tt[order];on=on[order];off=off[order]
    if np.any(np.diff(tt)<=0):raise ValueError('Duplicate/non-increasing extracted CAL timestamps')
    if tmin<tt[0] or tmax>tt[-1]:raise ValueError('Extracted CAL does not bracket selected science timestamps')
    onoff=(on-off)/20480*(20480+28672)/2
    smooth=onoff.copy() if n_t_cal==1 else uniform_filter1d(onoff,n_t_cal,axis=0,mode='nearest')
    # Retain the legacy interpolation-at-own-timestamps arithmetic as well.
    smooth=interp1d(tt,smooth,axis=0,fill_value='extrapolate')(tt)
    noise_f,noise=diode_reference(str(Path(tcal_path).resolve()),beam)
    diode=interp1d(noise_f,noise,axis=1,bounds_error=True)(cf)
    with np.errstate(divide='ignore',invalid='ignore'):
        temperature=(diode[None]/smooth*off).astype('float32').astype(float)
    pad=max(160,(n_t+1)//2+2)
    a=max(0,bracket(tt,tmin)-pad);b=min(len(tt),bracket(tt,tmax,True)+pad)
    temperature=temperature[a:b];time=tt[a:b]
    valid=np.isfinite(temperature)
    if n_t==1:filtered=np.where(valid,temperature,np.nan)
    elif valid.all():filtered=uniform_filter1d(temperature,n_t,axis=0,mode='nearest')
    else:
        support=uniform_filter1d(valid.astype(float),n_t,axis=0,mode='nearest')
        total=uniform_filter1d(np.where(valid,temperature,0),n_t,axis=0,mode='nearest')
        filtered=np.divide(total,support,out=np.full_like(total,np.nan),where=support>=.9)
    return dict(mjd=time,frequency_ghz=cf,raw=temperature,smoothed=filtered,
        info=dict(origin='Regenerated in memory from extracted CAL; original Tsys file unavailable',CAL_files=list(paths),
            raw_PSR_extraction=False,duty_factor=1.2,onoff_smoothing_rows=n_t_cal,onoff_control_applied=True,temperature_smoothing_rows=n_t,
            tcal_smoothing_channels=16,tcal_path=str(tcal_path),input_time_rows=len(tt),rows_read=len(time),frequency_channels_read=len(cf),
            CAL_time_range=[float(tt[0]),float(tt[-1])],unavailable_original_reference_comparison=True))


def read_reference(ctx,tmin,tmax,science_freq,cfg):
    n_t=smoothing_rows(cfg.get('n_t_apply',cfg.get('n_t',10)))
    n_t_cal=smoothing_rows(cfg.get('n_t_cal',10))
    if Path(ctx['reference']).is_file() and not (cfg.get('extracted_cal_root') and cfg.get('tcal_path')):
        pack=reference_data(ctx['reference'],tmin,tmax,science_freq,n_t=n_t)
        pack['info'].update(onoff_control_applied=False,onoff_smoothing_rows=cfg.get('reference_onoff_rows'),origin='Precomputed Tsys; ON-OFF smoothing embedded in input')
        return pack
    if not cfg.get('extracted_cal_root') or not cfg.get('tcal_path'):
        raise FileNotFoundError('Missing Tsys reference: '+ctx['reference']+'. Configure extracted_cal_root and tcal_path to regenerate from extracted CAL.')
    entries=cal_index(cfg['extracted_cal_root'],ctx['scan'],ctx['beam'])
    # > two half-filter contexts, plus whole-file boundary allowance.
    extra=0.
    context_rows=n_t+n_t_cal
    if context_rows>600:
        nearest=min(entries,key=lambda pair:abs(pair[0]-(tmin+tmax)/2))[1]
        times=extracted_band(nearest,float(science_freq[0]),float(science_freq[-1]))[0]
        extra=(context_rows-600)/2*float(np.median(np.diff(times)))
    paths=[name for mjd,name in entries if tmin-120/86400-extra<=mjd<=tmax+100/86400+extra]
    if not paths:raise ValueError('No extracted CAL files around selected time')
    return rebuild_reference(paths,tmin,tmax,science_freq,cfg['tcal_path'],ctx['beam'],n_t=n_t,n_t_cal=n_t_cal)


def reference_bounds(ctx,cfg):
    if Path(ctx['reference']).is_file():
        with fixed_table(ctx['reference'],extension=2,random_access=True) as (table,_):
            return float(table['MJD'][0]),float(table['MJD'][-1])
    entries=cal_index(cfg['extracted_cal_root'],ctx['scan'],ctx['beam'])
    with fits.open(entries[0][1],memmap=False) as h:lo=float(h[2].data['MJD'][0])
    with fits.open(entries[-1][1],memmap=False) as h:hi=float(h[2].data['MJD'][-1])
    return lo,hi
