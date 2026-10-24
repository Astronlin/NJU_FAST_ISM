"""Masked full-native-band broad sine/cosine fitting, with per-row amplitudes."""
import numpy as np
from .broad_sinusoid import fit_broad_sinusoid
from .rfi_baseline import frequency_mask,C

def fit_global_broad(data,freq_ghz,scales,rest,intervals,fit_bin=128,protect=250.,extra_exclusions=()):
    a=np.asarray(data,dtype=float);f=np.asarray(freq_ghz);scales=np.asarray(scales)
    velocity=(1-f[None,:]*scales[:,None]/rest)*C
    excluded=frequency_mask(f,intervals)|np.any(abs(velocity)<=protect,axis=0)
    for lo,hi in extra_exclusions:excluded|=np.any((velocity>=lo)&(velocity<=hi),axis=0)
    # A shared period per short chunk reduces nonlinear-search cost. Each
    # spectrum retains independent constant, cosine and sine amplitudes.
    _,_,report,_=fit_broad_sinusoid(a.mean(axis=0),f,excluded=excluded,fit_bin=fit_bin,grid_size=61)
    phase=2*np.pi*(f-f.mean())*1000/report['period_MHz']
    basis=np.column_stack([np.ones(len(f)),np.cos(phase),np.sin(phase)])
    n=len(f)//fit_bin;stop=n*fit_bin;design=basis[:stop].reshape(n,fit_bin,3).mean(axis=1)
    yy=a[:,:stop].reshape(len(a),n,fit_bin).mean(axis=2)
    allowed=(~excluded[:stop]).reshape(n,fit_bin).all(axis=1)
    good=np.broadcast_to(allowed,yy.shape)&np.isfinite(yy)
    w=good.astype(float)
    for _ in range(8):
        matrix=np.einsum('rn,ni,nj->rij',w,design,design)
        target=np.einsum('rn,rn,ni->ri',w,yy,design)
        coeff=np.linalg.solve(matrix,target[...,None])[...,0]
        residual=yy-coeff@design.T
        sigma=np.array([max(.05,1.4826*np.median(abs(z[g]-np.median(z[g])))) for z,g in zip(residual,good)])
        w=good/np.sqrt(1+(residual/sigma[:,None])**2)
    model=coeff@basis.T
    report.update(period_policy='One nonlinear period per <=32-spectrum chunk; independent constant/cosine/sine per native spectrum.',
        protected_policy='Union of per-spectrum +/-'+str(protect)+' km/s and all identified RFI; complete native fit bins only.',
        median_per_spectrum_anchor_RMS_K=float(np.median([np.sqrt(np.mean(z[g]**2)) for z,g in zip(residual,good)])))
    return model,coeff,report
