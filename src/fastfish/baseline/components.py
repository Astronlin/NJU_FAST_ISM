"""Ordered, masked spectral baseline components; compact per-polarization fits."""
import numpy as np
C=299792.458  # km/s; keep the Qt client free of server fitting dependencies


def normalize_components(value):
    if not isinstance(value,list) or len(value)>16:raise ValueError('Baseline components must be a list of at most 16 entries')
    result=[]
    for item in value:
        if not isinstance(item,dict):raise ValueError('Each baseline component must be an object')
        kind=item.get('function','polynomial')
        if kind not in ('polynomial','sincos'):raise ValueError('Function must be polynomial or sincos')
        window=item.get('range_kms')
        if window is not None:
            window=list(map(float,window))
            if len(window)!=2 or not np.isfinite(window).all() or window[0]>=window[1]:raise ValueError('Component range must be [minimum,maximum] km/s, or null for the full native band')
        spec=dict(function=kind,range_kms=window)
        if kind=='polynomial':
            degree=item.get('degree',1)
            if isinstance(degree,bool) or int(degree)!=degree or not 0<=degree<=20:raise ValueError('Polynomial degree must be an integer in 0..20')
            spec['degree']=int(degree)
        else:
            period=item.get('period_MHz')
            if period is not None and (not np.isfinite(period) or period<=0):raise ValueError('Sine/cosine period must be positive MHz, or null for automatic search')
            spec['period_MHz']=None if period is None else float(period)
            bounds=item.get('period_range_MHz')
            if bounds is not None:
                bounds=list(map(float,bounds))
                if len(bounds)!=2 or not np.isfinite(bounds).all() or not 0<bounds[0]<bounds[1]:raise ValueError('Invalid sine/cosine period search interval')
                if period is not None:raise ValueError('Choose a fixed period or a search interval, not both')
                spec['period_range_MHz']=bounds
        result.append(spec)
    return result


def default_components(cfg):
    if cfg.get('baseline_components') is not None:return normalize_components(cfg['baseline_components'])
    limit=float(cfg.get('local_fit_kms',1500))
    return [dict(function='sincos',range_kms=None,period_MHz=None),dict(function='polynomial',degree=1,range_kms=[-limit,limit])]


def _basis(f,v,spec,period,normalization):
    if spec['function']=='polynomial':
        x=(v-normalization[0])/normalization[1]
        return np.polynomial.chebyshev.chebvander(x,spec['degree'])
    phase=2*np.pi*(f-f.mean())*1000/period
    return np.broadcast_to(np.column_stack([np.ones(len(f)),np.cos(phase),np.sin(phase)]),(len(v),len(f),3))


def evaluate_record(record,f,scales,rest,rows=slice(None)):
    scales=np.asarray(scales);v=(1-np.asarray(f)[None,:]*scales[:,None]/rest)*C
    values=np.empty((len(scales),len(f),2))
    for p in range(2):
        coeff=np.asarray(record['coefficients'][p],dtype=float)[rows]
        basis=_basis(f,v,record['spec'],record['period_MHz'][p],record['normalization_kms'])
        values[:,:,p]=np.einsum('rnk,rk->rn',basis,coeff)
    window=record['spec']['range_kms']
    if window is not None:values=np.where(((v>=window[0])&(v<=window[1]))[:,:,None],values,0.)
    return values


def reconstruct_components(pol,f,scales,cfg,records,rows=slice(None),return_components=False):
    residual=np.array(pol,dtype=float,copy=True);models=[]
    for record in records:
        model=evaluate_record(record,f,scales,cfg['rest_freq_ghz'],rows);residual-=model
        if return_components:models.append(model)
    return residual,models


def fit_components(pol,f,scales,intervals,cfg):
    """Fit in list order. Whole 128-channel anchors; no H I/RFI fit leakage.

    Each range is BOTH its fitting and subtraction domain. No component is
    extrapolated beyond that range. Degrees use a stable Chebyshev basis.
    Nonlinear periods are shared only within a short chunk/polarization;
    every timestamp and AA/BB has separate linear coefficients.
    """
    from fastfish.rfi.detection import frequency_mask
    from fastfish.baseline.sinusoid import fit_broad_sinusoid
    specs=normalize_components(cfg.get('baseline_components',default_components(cfg))) if cfg.get('baseline',True) else []
    f=np.asarray(f);scales=np.asarray(scales);residual=np.array(pol,dtype=float,copy=True)
    v=(1-f[None,:]*scales[:,None]/cfg['rest_freq_ghz'])*C;fit_bin=int(cfg.get('baseline_fit_bin',128))
    n=len(f)//fit_bin;stop=n*fit_bin
    if n<1:raise ValueError('Not enough native channels for baseline anchors')
    rfi=frequency_mask(f[:stop],intervals).reshape(n,fit_bin).any(axis=1)
    records=[];reports=[]
    for ci,spec in enumerate(specs):
        window=spec['range_kms'];low,high=(float(v.min()),float(v.max())) if window is None else window
        normalization=[(low+high)/2,(high-low)/2]
        if normalization[1]<=0:raise ValueError('Empty baseline range')
        allowed=(abs(v[:,:stop])>float(cfg['protect_kms']))&~np.broadcast_to(np.repeat(rfi,fit_bin),v[:,:stop].shape)
        if window is not None:allowed&=(v[:,:stop]>=low)&(v[:,:stop]<=high)
        for lo,hi in cfg.get('baseline_exclusions_kms',[]):allowed&=~((v[:,:stop]>=lo)&(v[:,:stop]<=hi))
        whole=allowed.reshape(len(v),n,fit_bin).all(axis=2)
        record=dict(spec=spec,normalization_kms=normalization,period_MHz=[],coefficients=[]);details=[]
        for p in range(2):
            period=spec.get('period_MHz');search={}
            if spec['function']=='sincos' and period is None:
                excluded=np.ones(len(f),dtype=bool);excluded[:stop]=~np.repeat(whole.all(axis=0),fit_bin)
                _,_,search,_=fit_broad_sinusoid(np.mean(residual[:,:,p],axis=0),f,excluded=excluded,fit_bin=fit_bin,grid_size=61,period_range_mhz=spec.get('period_range_MHz'))
                period=search['period_MHz']
            if spec['function']=='sincos' and period<4*fit_bin*abs(np.median(np.diff(f)))*1000:raise ValueError('Sine/cosine period is undersampled by the fitting anchors')
            basis=_basis(f,v,spec,period,normalization);k=basis.shape[-1];design=basis[:,:stop].reshape(len(v),n,fit_bin,k).mean(axis=2)
            anchors=residual[:,:stop,p].reshape(len(v),n,fit_bin).mean(axis=2);good=whole&np.isfinite(anchors)
            required=max(8,3*k)
            if (good.sum(axis=1)<required).any():raise ValueError('Component %d (%s): need at least %d clean anchors per spectrum'%(ci+1,spec['function'],required))
            vv=v[:,:stop].reshape(len(v),n,fit_bin).mean(axis=2)
            if low < -cfg['protect_kms'] and high > cfg['protect_kms']:
                if ((good&(vv< -cfg['protect_kms'])).sum(axis=1)<8).any() or ((good&(vv>cfg['protect_kms'])).sum(axis=1)<8).any():raise ValueError('Component %d: need eight clean anchors on each side of protected H I'%(ci+1))
            target=np.where(good,anchors,0.);w=good.astype(float)
            for _ in range(8):
                matrix=np.einsum('rn,rni,rnj->rij',w,design,design);rhs=np.einsum('rn,rn,rni->ri',w,target,design)
                if not np.isfinite(matrix).all() or np.any(np.linalg.cond(matrix)>1e12):raise ValueError('Component %d fit is ill-conditioned; lower the degree or revise its range'%(ci+1))
                coeff=np.linalg.solve(matrix,rhs[...,None])[...,0];rr=target-np.einsum('rnk,rk->rn',design,coeff)
                scatter=np.array([max(.05,1.4826*np.median(abs(a[g]-np.median(a[g])))) for a,g in zip(rr,good)])
                w=good/np.sqrt(1+(rr/scatter[:,None])**2)
            record['coefficients'].append(coeff.tolist());record['period_MHz'].append(None if period is None else float(period))
            details.append(dict(polarization=['AA','BB'][p],minimum_anchors=int(good.sum(axis=1).min()),median_anchor_RMS_K=float(np.median([np.sqrt(np.mean(a[g]**2)) for a,g in zip(rr,good)])),period_search=search if search else None))
        residual-=evaluate_record(record,f,scales,cfg['rest_freq_ghz']);records.append(record);reports.append(dict(component=ci+1,spec=spec,polarizations=details))
    return records,dict(baseline=bool(cfg.get('baseline',True)),baseline_mode='ordered_components',components=reports,at_search_boundary=any(p['period_search'] and p['period_search'].get('at_search_boundary',False) for r in reports for p in r['polarizations']))
