"""Read-only, bounded preview engine shared by interactive and batch stages."""
from pathlib import Path
from collections import OrderedDict
import hashlib
import json
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
from fastfish.calibration.reference import reference_data,smoothing_rows
from fastfish.calibration.extracted import read_reference
from fastfish.io.fits_stream import fixed_table
from fastfish.processing.spectral import (frequency_gain, affine_overlap_rebin,
    gaussian_to_fwhm, affine_linear, channel_edges)
from fastfish.rfi.detection import C, track_intervals, frequency_mask, fit_linear_anchors
from fastfish.baseline.broad import fit_global_broad
from fastfish.baseline.components import fit_components,reconstruct_components
from fastfish.baseline.windows import normalize_windows
from fastfish.rfi.catalog import load_catalog, profile_priors, catalog_fingerprint, catalog_path


def calibrate(raw, f, reference, cf, scales, priors, cfg, compact=False):
    anchors=raw.mean(axis=(0,2)).reshape(-1,128).mean(axis=1)
    intervals,tracks=track_intervals(anchors,f.reshape(-1,128).mean(axis=1),priors,float(np.median(scales)),cfg['rest_freq_ghz'],protect=float(cfg.get('protect_kms',250.)),catalog=cfg.get('rfi_catalog','packaged'))
    rest=cfg['rest_freq_ghz']
    excluded=[((1-700/C)*rest,(1+700/C)*rest)]+[(a/1000,b/1000) for a,b in intervals]
    gain,details=frequency_gain(raw,f,reference,cf,excluded_ranges=excluded,
                               degree=cfg['gain_degree'],return_coefficients=True,evaluate=not compact)
    if compact:
        return dict(coefficients=details['coefficients'],reference_edges_ghz=details['reference_edges_ghz'],intervals=intervals,tracks=tracks)
    pol=raw/gain
    return dict(pol=pol,average=pol.mean(axis=2),coefficients=details['coefficients'],
                reference_edges_ghz=details['reference_edges_ghz'],intervals=intervals,tracks=tracks)


def reconstruct_calibrated(raw, f, coefficients, reference_edges):
    x=(f-(f[0]+f[-1])/2)/((f[-1]-f[0])/2)
    gain=np.exp(np.polynomial.polynomial.polyval(x,coefficients.transpose(2,0,1))).transpose(0,2,1)
    supported=(f>=reference_edges[0])&(f<=reference_edges[1])
    gain[:,~supported,:]=np.nan
    pol=raw/gain
    return pol


def local_linear_correction(coefficients, velocity, limit):
    """Apply a fitted local correction only within its closed LSRK window."""
    coefficients=np.asarray(coefficients)
    velocity=np.asarray(velocity)
    correction=coefficients[:,0,None]+coefficients[:,1,None]*velocity/limit
    return np.where(np.abs(velocity)<=limit,correction,0.)


def fit_baseline(average, f, scales, intervals, cfg, sensitivity=False, compact=False):
    extra=normalize_windows(cfg.get('baseline_exclusions_kms',[]))
    protect=float(cfg['protect_kms']);limit=float(cfg['local_fit_kms'])
    if not 0<=protect<limit:
        raise ValueError('Local fit half-width must exceed the protected half-width')
    velocity=(1-f[None,:]*scales[:,None]/cfg['rest_freq_ghz'])*C
    if cfg.get('baseline',True):
        broad,gcoeff,report=fit_global_broad(average,f,scales,cfg['rest_freq_ghz'],intervals,protect=protect,extra_exclusions=extra)
        coeff,info=fit_linear_anchors(average-broad,f,scales,cfg['rest_freq_ghz'],intervals,fit_limit=limit,protect=protect,extra_exclusions=extra)
        local=None if compact else local_linear_correction(coeff,velocity,limit)
        alt=None
        if sensitivity:
            alt_limit=float(cfg['alternative_fit_kms'])
            alt,_=fit_linear_anchors(average-broad,f,scales,cfg['rest_freq_ghz'],intervals,fit_limit=alt_limit,protect=protect,extra_exclusions=extra)
    else:
        broad=np.zeros_like(average);local=np.zeros_like(average)
        coeff=np.zeros((len(average),2));gcoeff=np.zeros((len(average),3));alt=coeff.copy() if sensitivity else None
        report=dict(baseline=False,period_MHz=1.);info={}
    return dict(broad=None if compact else broad,local=local,residual=None if compact else average-broad-local,
                global_coeff=gcoeff,local_coeff=coeff,alternative_coeff=alt,report=report,info=info,local_fit_kms=limit)


def local_polarization_mean(coefficients, velocity, limit):
    return sum(local_linear_correction(coefficients[:,p,:],velocity,limit) for p in range(2))/2


def fit_polarization_baselines(pol,f,scales,intervals,cfg,sensitivity=False,compact=False):
    """Fit each auto-polarization independently, then average corrected spectra."""
    if pol.ndim!=3 or pol.shape[2]!=2:
        raise ValueError('Expected calibrated AA/BB with shape (time, frequency, 2)')
    fits=[fit_baseline(pol[:,:,p],f,scales,intervals,cfg,sensitivity,compact=compact) for p in range(2)]
    result={key:np.stack([fit[key] for fit in fits],axis=2) for key in (() if compact else ('broad','local','residual'))}
    result={**{key:values.mean(axis=2) for key,values in result.items()},
            **{key+'_pol':values for key,values in result.items()}}
    for key in ('global_coeff','local_coeff','alternative_coeff'):
        result[key]=np.stack([fit[key] for fit in fits],axis=1) if fits[0][key] is not None else None
    result.update(report=dict(baseline=bool(cfg.get('baseline',True)),baseline_mode='per_polarization',
        polarizations={label:fit['report'] for label,fit in zip(('AA','BB'),fits)},
        at_search_boundary=any(fit['report'].get('at_search_boundary',False) for fit in fits)),
        info={label:fit['info'] for label,fit in zip(('AA','BB'),fits)},local_fit_kms=float(cfg['local_fit_kms']))
    return result


def reconstruct_polarization_baselines(pol,f,scales,global_coeff,local_coeff,periods,cfg,alternative_coeff=None):
    """Reconstruct compact AA/BB fits in bounded memory; return their corrected mean."""
    if global_coeff.shape!=(len(pol),2,3) or local_coeff.shape!=(len(pol),2,2) or np.shape(periods)!=(2,):
        raise ValueError('Baseline records must contain independent AA/BB coefficients and periods')
    velocity=(1-f[None,:]*scales[:,None]/cfg['rest_freq_ghz'])*C
    residual=np.zeros(pol.shape[:2]);delta=np.zeros_like(residual)
    for p in range(2):
        phase=2*np.pi*(f-f.mean())*1000/periods[p]
        gc=global_coeff[:,p,:]
        broad=gc[:,0,None]+gc[:,1,None]*np.cos(phase)+gc[:,2,None]*np.sin(phase)
        local=local_linear_correction(local_coeff[:,p,:],velocity,cfg['local_fit_kms'])
        residual+=(pol[:,:,p]-broad-local)/2
        if alternative_coeff is not None:
            delta+=(local-local_linear_correction(alternative_coeff[:,p,:],velocity,cfg['alternative_fit_kms']))/2
    return residual,delta


def output_grid(cfg):
    lo,hi=map(float,cfg['velocity_range_kms'])
    if not np.isfinite([lo,hi]).all() or lo>=hi:
        raise ValueError('Product velocity minimum must be below maximum')
    width=float(cfg['n_b'] if cfg['spectral_mode']=='rebin' else cfg['n_g'])
    if width<1 or cfg['n_g']<cfg['native_fwhm_chan']:
        raise ValueError('Invalid rebin/Gaussian width')
    step=cfg['native_df_hz']*width/1e9;rest=cfg['rest_freq_ghz']
    lower=(1-hi/C)*rest;upper=(1-lo/C)*rest
    origin=cfg.get('output_frequency_origin_ghz',lower)
    first=origin+np.ceil((lower-origin)/step-1e-9)*step
    q=np.arange(first,upper,step)
    if len(q)<2:
        raise ValueError('Product range contains fewer than two output channels')
    return q


def resample(residual, f, scales, intervals, cfg, q=None):
    q=output_grid(cfg) if q is None else q
    mark=frequency_mask(f,intervals,half_channel=cfg['native_df_hz']/2e6).astype(float)
    marks=np.broadcast_to(mark,residual.shape)
    if cfg['spectral_mode']=='rebin':
        values,support=affine_overlap_rebin(residual,f,q,scales,axis=1,min_coverage=1.)
        fraction=affine_overlap_rebin(marks,f,q,scales,axis=1,min_coverage=1.)[0]
    else:
        values=affine_linear(gaussian_to_fwhm(residual,cfg['n_g'],cfg['native_fwhm_chan']),f,q,scales,axis=1)
        fraction=affine_linear(gaussian_to_fwhm(marks,cfg['n_g'],cfg['native_fwhm_chan']),f,q,scales,axis=1)
    if not np.isfinite(values).all() or not np.isfinite(fraction).all():
        raise ValueError('Output range or inputs lack complete finite spectral support')
    return values,np.clip(fraction,0,1)


def _bytes(value):
    if isinstance(value,np.ndarray):return value.nbytes+128
    if isinstance(value,dict):return sys.getsizeof(value)+sum(_bytes(k)+_bytes(v) for k,v in value.items())
    if isinstance(value,(list,tuple)):return sys.getsizeof(value)+sum(_bytes(v) for v in value)
    return sys.getsizeof(value)


class PreviewDataset:
    """Reads inputs/selection only. Never configures a run or writes a product."""
    def __init__(self, config_path, selection_path, memory_mib=512):
        self.config_path=Path(config_path).resolve();self.selection_path=Path(selection_path).resolve()
        self.cfg=json.loads(self.config_path.read_text());self.cfg.setdefault('n_t_apply',self.cfg.get('n_t',10));self.cfg.setdefault('n_t_cal',10);self.selection=json.loads(self.selection_path.read_text())
        catalog=load_catalog(self.cfg.get('rfi_catalog','packaged'),self.config_path.parent)
        self.cfg['rfi_catalog']=str(catalog_path(self.cfg.get('rfi_catalog','packaged'),self.config_path.parent))
        self.catalog_sha256=catalog_fingerprint(self.cfg.get('rfi_catalog','packaged'),self.config_path.parent)
        self.priors=profile_priors(catalog,self.cfg['rfi_profile'],self.cfg['scans']) if self.cfg.get('rfi_profile') else {scan:[] for scan in self.cfg['scans']}
        self.files=[];coordinates=[];self.offsets=[0]
        for item in self.selection['items']:
            for ctx in item['files']:
                with np.load(ctx['coordinates']) as p:c={k:p[k].copy() for k in p.files}
                self.files.append(dict(ctx,scan=item['scan'],beam=item['beam'],reference=item['reference'],coords=c))
                coordinates.append(c);self.offsets.append(self.offsets[-1]+len(c['ids']))
        self.offsets=np.asarray(self.offsets)
        self.ra=np.concatenate([c['ra'] for c in coordinates]);self.dec=np.concatenate([c['dec'] for c in coordinates])
        self.groups=np.concatenate([np.full(len(x['coords']['ids']),self.cfg['scans'].index(x['scan']),dtype='uint8') for x in self.files])
        self.tree_dec=float(np.median(self.dec))
        self.tree=cKDTree(np.column_stack([self.ra*np.cos(np.deg2rad(self.tree_dec)),self.dec]))
        self.cache=OrderedDict();self.cache_bytes=0;self.memory_limit=int(memory_mib*1024**2)
        self.input_signature=hashlib.sha256(self.config_path.read_bytes()+self.selection_path.read_bytes()).hexdigest()

    def address(self,index):
        file_id=int(np.searchsorted(self.offsets,index,side='right')-1)
        local=int(index-self.offsets[file_id]);size=self.cfg['chunk_rows']
        return file_id,local//size,local%size

    def closest(self,ra,dec):
        return int(self.tree.query([ra*np.cos(np.deg2rad(self.tree_dec)),dec])[1])

    def rectangle(self,x0,x1,y0,y1):
        return np.flatnonzero((self.ra>=min(x0,x1))&(self.ra<=max(x0,x1))&(self.dec>=min(y0,y1))&(self.dec<=max(y0,y1)))

    def polygon(self,vertices):
        points=np.asarray(vertices,dtype=float)
        if points.ndim!=2 or points.shape[1]!=2 or len(points)<3 or not np.isfinite(points).all():
            raise ValueError('Polygon requires at least three finite RA/Dec vertices')
        ids=self.rectangle(points[:,0].min(),points[:,0].max(),points[:,1].min(),points[:,1].max())
        x=self.ra[ids];y=self.dec[ids];inside=np.zeros(len(ids),dtype=bool);edge=inside.copy()
        for (x0,y0),(x1,y1) in zip(points,np.roll(points,-1,axis=0)):
            dx=x1-x0;dy=y1-y0
            if dx==0 and dy==0:continue
            cross=(x-x0)*dy-(y-y0)*dx
            edge|=(np.abs(cross)<=1e-12*max(1,abs(dx),abs(dy)))&(x>=min(x0,x1))&(x<=max(x0,x1))&(y>=min(y0,y1))&(y<=max(y0,y1))
            if dy!=0:inside^=((y0>y)!=(y1>y))&(x<x0+(y-y0)*dx/dy)
        return ids[inside|edge]

    def _cached(self,key,compute):
        if key in self.cache:
            value=self.cache.pop(key);self.cache[key]=value;return value
        value=compute();size=_bytes(value)
        # An uncachable large result must not evict reusable small entries.
        if size>self.memory_limit:return value
        while self.cache and self.cache_bytes+size>self.memory_limit:
            _,old=self.cache.popitem(last=False);self.cache_bytes-=_bytes(old)
        if size<=self.memory_limit:self.cache[key]=value;self.cache_bytes+=size
        return value

    def reference(self,file_id,chunk_id,params=None):
        n_t=smoothing_rows((params or self.cfg).get('n_t_apply',10));n_t_cal=smoothing_rows((params or self.cfg).get('n_t_cal',10))
        ctx=self.files[file_id];size=self.cfg['chunk_rows'];sl=slice(chunk_id*size,(chunk_id+1)*size)
        t=ctx['coords']['mjd'][sl];f=self.frequency(file_id)
        return self._cached(('reference',file_id,chunk_id,n_t,n_t_cal),lambda:read_reference(ctx,float(t.min()),float(t.max()),f,dict(self.cfg,n_t_apply=n_t,n_t_cal=n_t_cal)))

    def frequency(self,file_id):
        e=self.files[file_id]['entry'];key=('frequency',e['freq_first_ghz'],e['channels'],e['native_df_hz'])
        return self._cached(key,lambda:e['freq_first_ghz']+np.arange(e['channels'])*e['native_df_hz']/1e9)

    def raw(self,file_id,chunk_id):
        """SCI counts are independent of the two reference smoothing controls."""
        def compute():
            ctx=self.files[file_id];size=self.cfg['chunk_rows'];sl=slice(chunk_id*size,(chunk_id+1)*size)
            with fixed_table(ctx['entry']['path'],random_access=True) as (tab,_):
                return tab['DATA'][ctx['coords']['ids'][sl],:,:2].astype(float)
        return self._cached(('raw',file_id,chunk_id),compute)

    def reference_display(self,file_id,chunk_id,params):
        """Diagnostic band summaries; filtering still acts on each channel."""
        key=('reference_display',file_id,chunk_id,smoothing_rows(params.get('n_t_apply',10)),smoothing_rows(params.get('n_t_cal',10)))
        def compute():
            pack=self.reference(file_id,chunk_id,params)
            return dict(reference_time=pack['mjd'],reference_raw=np.nanmedian(pack['raw'],axis=2),reference_smooth=np.nanmedian(pack['smoothed'],axis=2))
        return self._cached(key,compute)

    def scales(self,file_id):
        """Exact per-spectrum Doppler factors, deferred until science processing."""
        c=self.files[file_id]['coords']
        if 'frequency_scale' in c:return c['frequency_scale']
        def compute():
            from fastfish.pointing.doppler import frequency_scales
            from fastfish.pointing.time_tables import initialize_earth_orientation
            initialize_earth_orientation()
            return frequency_scales(c['mjd'],c['ra'],c['dec'],mode='spectrum')
        return self._cached(('frequency_scales',file_id),compute)

    def calibrated(self,file_id,chunk_id,params=None):
        n_t=smoothing_rows((params or self.cfg).get('n_t_apply',10));n_t_cal=smoothing_rows((params or self.cfg).get('n_t_cal',10))
        def compute():
            ctx=self.files[file_id];size=self.cfg['chunk_rows'];sl=slice(chunk_id*size,(chunk_id+1)*size)
            f=self.frequency(file_id);c=ctx['coords'];pack=self.reference(file_id,chunk_id,dict(n_t_apply=n_t,n_t_cal=n_t_cal))
            ref=interp1d(pack['mjd'],pack['smoothed'],axis=0,bounds_error=True)(c['mjd'][sl])
            raw=self.raw(file_id,chunk_id)
            scales=self.scales(file_id)[sl]
            result=calibrate(raw,f,ref,pack['frequency_ghz'],scales,self.priors[ctx['scan']],self.cfg)
            result.update(raw_pol=raw,frequency_ghz=f,scales=scales,mjd=c['mjd'][sl])
            return result
        return self._cached(('calibrated',file_id,chunk_id,n_t,n_t_cal),compute)

    def baseline_key(self,file_id,chunk_id,params):
        return ('baseline',file_id,chunk_id,smoothing_rows(params.get('n_t_apply',10)),smoothing_rows(params.get('n_t_cal',10)),bool(params['baseline']),float(params['protect_kms']),float(params['local_fit_kms']),tuple(map(tuple,normalize_windows(params.get('baseline_exclusions_kms',[])))),json.dumps(params.get('baseline_components'),sort_keys=True))

    def baseline(self,file_id,chunk_id,params,rows=None,cal=None):
        key=self.baseline_key(file_id,chunk_id,params)
        cal=self.calibrated(file_id,chunk_id,params) if cal is None else cal
        if params.get('baseline_components') is not None:
            cfg=dict(self.cfg,**params)
            def component_compute():
                records,report=fit_components(cal['pol'],cal['frequency_ghz'],cal['scales'],cal['intervals'],cfg)
                return dict(component_records=records,report=report)
            compact=self._cached(key,component_compute);rows=slice(None) if rows is None else rows
            pol=cal['pol'][rows];residual,models=reconstruct_components(pol,cal['frequency_ghz'],cal['scales'][rows],cfg,compact['component_records'],rows,return_components=True)
            first=models[0] if models else np.zeros_like(pol);remaining=sum(models[1:],np.zeros_like(pol))
            return dict(report=compact['report'],component_records=compact['component_records'],component_models=models,
                broad_pol=first,local_pol=remaining,residual_pol=residual,broad=first.mean(axis=2),local=remaining.mean(axis=2),residual=residual.mean(axis=2))
        def compute():
            fitted=fit_polarization_baselines(cal['pol'],cal['frequency_ghz'],cal['scales'],cal['intervals'],dict(self.cfg,**params))
            # Retain coefficients/support diagnostics, not six native model cubes.
            return {k:fitted[k] for k in ('global_coeff','local_coeff','alternative_coeff','report','info','local_fit_kms')}
        compact=self._cached(key,compute)
        rows=slice(None) if rows is None else rows
        pol=cal['pol'][rows];f=cal['frequency_ghz'];scales=cal['scales'][rows]
        velocity=(1-f[None,:]*scales[:,None]/self.cfg['rest_freq_ghz'])*C
        gc=compact['global_coeff'][rows];lc=compact['local_coeff'][rows]
        broad=[];local=[]
        for p,label in enumerate(('AA','BB')):
            period=compact['report']['polarizations'][label]['period_MHz'];phase=2*np.pi*(f-f.mean())*1000/period
            broad.append(gc[:,p,0,None]+gc[:,p,1,None]*np.cos(phase)+gc[:,p,2,None]*np.sin(phase))
            local.append(local_linear_correction(lc[:,p,:],velocity,compact['local_fit_kms']))
        broad=np.stack(broad,axis=2);local=np.stack(local,axis=2);residual=pol-broad-local
        return dict(compact,global_coeff=gc,local_coeff=lc,broad_pol=broad,local_pol=local,residual_pol=residual,broad=broad.mean(axis=2),local=local.mean(axis=2),residual=residual.mean(axis=2))

    def point(self,index,params,selected_only=False):
        file_id,chunk_id,row=self.address(index);cal=self.calibrated(file_id,chunk_id,params)
        count=len(cal['mjd']);chunk_row=row
        result=self.baseline(file_id,chunk_id,params,slice(row,row+1) if selected_only else None,cal=cal)
        if selected_only:
            cal=dict(cal,**{key:cal[key][row:row+1] for key in ('raw_pol','pol','average','scales','mjd')});row=0
        return dict(index=int(index),file_id=file_id,chunk_id=chunk_id,row=row,chunk_row=chunk_row,chunk_size=count,fit_key=self.baseline_key(file_id,chunk_id,params),cal=cal,fit=result,baseline_exclusions_kms=normalize_windows(params.get('baseline_exclusions_kms',[])),reference=self.reference(file_id,chunk_id,params))

    def region_spectrum(self,indices,params,max_spectra=128,cancel=None,progress=None):
        ids=np.asarray(indices,dtype=int)
        if not len(ids):raise ValueError('Select a region containing spectra')
        total=len(ids)
        preferred=None
        if max_spectra and total>max_spectra:
            files=np.searchsorted(self.offsets,ids,side='right')-1
            chunks=(ids-self.offsets[files])//self.cfg['chunk_rows']
            keys=files.astype('int64')*1000000+chunks
            unique=np.unique(keys);maximum=max(1,int(np.ceil(max_spectra/self.cfg['chunk_rows'])))
            chosen=unique[np.unique(np.linspace(0,len(unique)-1,min(maximum,len(unique))).astype(int))].tolist()
            for key in reversed(self.cache):
                if key[0]=='calibrated':
                    candidate=key[1]*1000000+key[2]
                    if candidate in unique:
                        preferred=(key[1],key[2]);chosen=[candidate]+[x for x in chosen if x!=candidate][:maximum-1];break
            ids=ids[np.isin(keys,chosen)]
        jobs={}
        for index in ids:
            fi,ch,row=self.address(int(index));jobs.setdefault((fi,ch),[]).append(row)
        if preferred in jobs:jobs={preferred:jobs[preferred],**{k:v for k,v in jobs.items() if k!=preferred}}
        cfg=dict(self.cfg,**params);q=output_grid(cfg);summed=np.zeros(len(q));marked=np.zeros(len(q));count=0
        for (fi,ch),rows in jobs.items():
            if cancel is not None and cancel():return None
            cal=self.calibrated(fi,ch,params);fit=self.baseline(fi,ch,params,rows=rows,cal=cal)
            values,fraction=resample(fit['residual'],cal['frequency_ghz'],cal['scales'][rows],cal['intervals'],cfg,q)
            summed+=values.sum(axis=0);marked+=fraction.sum(axis=0);count+=len(rows)
            if progress is not None:
                progress(dict(velocity=(1-q/cfg['rest_freq_ghz'])*C,mean=summed/count,rfi_fraction=marked/count,used=count,total=total,target=len(ids)))
        return dict(velocity=(1-q/cfg['rest_freq_ghz'])*C,mean=summed/count,rfi_fraction=marked/count,used=count,total=total)
