"""Profile a broad sinusoid on the full supplied band, masking spectral lines."""
import numpy as np
from scipy.optimize import minimize_scalar


def fit_broad_sinusoid(data,freq_ghz,excluded=None,fit_bin=128,period_range_mhz=None,grid_size=401):
    """Return full-band model, coefficients and a period-profile diagnostic.

    Model: constant + a*cos(2*pi*(nu-centre)/P)+b*sin(...).
    Line exclusions preserve the full native coordinate and reject whole fit
    anchors crossing a boundary. Fit-only means use exact averaged sinusoidal
    basis functions. Soft-L1 iteratively reweighted fits and one fixed robust
    scale score the SAME initial anchors for every candidate period.

    Default broad-period range = one-quarter to four times the observed band.
    That range is a modelling choice, not a measured instrumental period.
    The near-minimum profile range is a sensitivity diagnostic, not a CI.
    """
    y=np.asarray(data,dtype=float);f=np.asarray(freq_ghz,dtype=float)
    if y.ndim!=1 or f.shape!=y.shape or len(f)<16 or not np.isfinite(f).all() or np.any(np.diff(f)<=0):raise ValueError('Expected one spectrum and increasing frequencies')
    if int(fit_bin)!=fit_bin or fit_bin<1 or int(grid_size)!=grid_size or grid_size<21:raise ValueError('Invalid bin width or search size')
    fit_bin=int(fit_bin);grid_size=int(grid_size);df=float(np.median(np.diff(f)))
    if not np.allclose(np.diff(f),df,rtol=1e-6,atol=abs(df)*1e-7):raise ValueError('The exact anchor response assumes uniform native spacing')
    nu=(f-(f[0]+f[-1])/2)*1000;span=(f[-1]-f[0])*1000
    if period_range_mhz is None:period_range_mhz=(span/4,span*4)
    pmin,pmax=map(float,period_range_mhz)
    if not np.isfinite([pmin,pmax]).all() or not 0<pmin<pmax or pmin<4*fit_bin*df*1000:raise ValueError('Invalid or undersampled broad-period range')
    allowed=np.isfinite(y)
    if excluded is not None:
        excluded=np.asarray(excluded,dtype=bool)
        if excluded.shape!=y.shape:raise ValueError('Exclusion mask must match spectrum')
        allowed&=~excluded
    n=len(f)//fit_bin;stop=n*fit_bin
    good=allowed[:stop].reshape(n,fit_bin).all(axis=1)
    if good.sum()<20:raise ValueError('Too few unmasked full-band anchors')
    centres=nu[:stop].reshape(n,fit_bin).mean(axis=1)
    anchor=y[:stop].reshape(n,fit_bin).mean(axis=1)
    xx=centres[good];yy=anchor[good]
    adjacent=good[:-1]&good[1:]
    differences=np.diff(anchor)[adjacent]
    if len(differences)<5:raise ValueError('Too few adjacent anchors to estimate robust scatter')
    sigma=max(1e-8,float(1.4826*np.median(np.abs(differences-np.median(differences)))/np.sqrt(2)))
    def design(period):
        response=np.sinc(fit_bin*df*1000/period)/np.sinc(df*1000/period)
        phase=2*np.pi*xx/period
        return np.column_stack([np.ones(len(xx)),response*np.cos(phase),response*np.sin(phase)])
    def solve(cycles):
        period=span/cycles;basis=design(period);coeff=np.linalg.lstsq(basis,yy,rcond=None)[0]
        for _ in range(20):
            residual=yy-basis@coeff;weight=1/np.sqrt(1+(residual/sigma)**2);root=np.sqrt(weight)
            update=np.linalg.lstsq(basis*root[:,None],yy*root,rcond=None)[0]
            if np.max(np.abs(update-coeff))<1e-10:coeff=update;break
            coeff=update
        residual=yy-basis@coeff
        cost=float(np.sum(2*sigma*sigma*(np.sqrt(1+(residual/sigma)**2)-1)))
        return cost,coeff,residual,float(np.linalg.cond(basis))
    cycles=np.linspace(span/pmax,span/pmin,grid_size)
    # All trial periods share the same clean anchors. Batch the coarse search;
    # refinement still uses the original least-squares solver independently.
    periods=span/cycles
    response=np.sinc(fit_bin*df*1000/periods)/np.sinc(df*1000/periods)
    phase=2*np.pi*xx[None,:]/periods[:,None]
    designs=np.stack([np.ones_like(phase),response[:,None]*np.cos(phase),response[:,None]*np.sin(phase)],axis=2)
    weights=np.ones(phase.shape);trial=np.zeros((len(cycles),3));active=np.ones(len(cycles),dtype=bool)
    # One unweighted solve followed by up to twenty IRLS updates, matching solve().
    for iteration in range(21):
        matrix=np.einsum('pn,pni,pnj->pij',weights,designs,designs)
        target=np.einsum('pn,pni,n->pi',weights,designs,yy)
        updated=np.linalg.solve(matrix,target[...,None])[...,0]
        converged=(np.max(np.abs(updated-trial),axis=1)<1e-10) if iteration else np.zeros(len(cycles),dtype=bool)
        trial[active]=updated[active]
        active&=~converged
        residuals=yy[None,:]-np.einsum('pni,pi->pn',designs,trial)
        weights=1/np.sqrt(1+(residuals/sigma)**2)
        if not active.any():break
    cost=np.sum(2*sigma*sigma*(np.sqrt(1+(residuals/sigma)**2)-1),axis=1)
    best=int(np.argmin(cost))
    candidates=[(float(cycles[best]),solve(cycles[best]))]
    # Refine every discrete local minimum, avoiding a single basin assumption.
    minima=np.flatnonzero((cost[1:-1]<=cost[:-2])&(cost[1:-1]<=cost[2:]))+1
    for j in minima:
        opt=minimize_scalar(lambda x:solve(x)[0],bounds=(cycles[j-1],cycles[j+1]),method='bounded',options={'xatol':1e-10})
        candidates.append((float(opt.x),solve(opt.x)))
    selected,record=min(candidates,key=lambda x:x[1][0]);best_cost,coeff,residual,condition=record
    if condition>1e6:raise ValueError('Broad sinusoid is poorly conditioned')
    period=span/selected;phase=2*np.pi*nu/period;model=coeff[0]+coeff[1]*np.cos(phase)+coeff[2]*np.sin(phase)
    near=cost<=best_cost+max(.01*abs(best_cost),sigma*sigma)
    report={'method':'Full-band constant plus broad sine/cosine; soft-L1 profile fit',
        'period_MHz':float(period),'amplitude_K':float(np.hypot(coeff[1],coeff[2])),
        'constant_K':float(coeff[0]),'cosine_K':float(coeff[1]),'sine_K':float(coeff[2]),
        'phase_origin_GHz':float((f[0]+f[-1])/2),'native_bandwidth_MHz':float(span),
        'period_search_range_MHz':[pmin,pmax],'cycles_across_band':float(selected),
        'at_search_boundary':bool(abs(selected-cycles[0])<1e-5 or abs(selected-cycles[-1])<1e-5),
        'fit_anchors':int(good.sum()),'total_anchors':int(n),'fit_bin_native_channels':fit_bin,
        'robust_anchor_scale_K':sigma,'design_condition':condition,'profile_cost':best_cost,
        'anchors_beyond_4_sigma':int((np.abs(residual)>4*sigma).sum()),
        'near_minimum_profile_period_range_MHz':[float(np.min(span/cycles[near])),float(np.max(span/cycles[near]))] if near.any() else [float(period),float(period)],
        'profile_range_definition':'Cost within max(1 percent of optimum cost, anchor scale squared); not a statistical confidence interval.',
        'full_band_mean_anchor_residual_K':float(np.mean(residual)),
        'full_band_median_anchor_residual_K':float(np.median(residual))}
    diagnostic={'anchor_freq_ghz':centres/1000+(f[0]+f[-1])/2,'anchor_temperature_K':anchor,
        'anchor_used':good,'trial_periods_MHz':span/cycles,'trial_cost':cost}
    return model,coeff,report,diagnostic
