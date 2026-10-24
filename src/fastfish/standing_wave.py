"""Conditional, H I-protected short-period residual correction.

A held-out frequency-block check decides whether a chunk contains a coherent
~1.09 MHz component. Each retained spectrum/polarization then has independent
robust linear/sine/cosine coefficients. No fit uses protected H I.
"""
import numpy as np
from .rfi_baseline import C,frequency_mask


def evaluate_wave(coefficients,period,freq_ghz,velocity,limit):
    phase=2*np.pi*(np.asarray(freq_ghz)-np.mean(freq_ghz))*1000/period
    values=coefficients[:,0,None]*np.cos(phase)+coefficients[:,1,None]*np.sin(phase)
    return np.where(abs(velocity)<=limit,values,0.)


def fit_residual_wave(residual,freq_ghz,scales,rest,intervals,limit=1500.,protect=250.,fit_bin=128):
    """Return additional local coefficients and conditional ripple coefficients.

    Selection uses a quadratic nuisance trend (not subtracted) and alternating
    blocks of eight fit anchors. Require >=10% held-out RMS improvement and
    amplitude >=0.03 K. These explicit selection choices require validation
    for other observing setups. A failure to detect returns exactly zero.
    """
    y=np.asarray(residual);f=np.asarray(freq_ghz);scales=np.asarray(scales)
    n=len(f)//fit_bin;stop=n*fit_bin
    ff=f[:stop].reshape(n,fit_bin);anchors=y[:,:stop].reshape(len(y),n,fit_bin).mean(axis=2)
    v=(1-ff.mean(axis=1)[None,:]*scales[:,None]/rest)*C
    vl=(1-ff[:,0][None,:]*scales[:,None]/rest)*C;vr=(1-ff[:,-1][None,:]*scales[:,None]/rest)*C
    good=(vl<=limit)&(vr>=-limit)&((vr>protect)|(vl<-protect))&~frequency_mask(f[:stop],intervals).reshape(n,fit_bin).any(axis=1)[None,:]&np.isfinite(anchors)
    common=good.all(axis=0);index=np.arange(n);train=common&((index//8)%2==0);test=common&~train
    local=np.zeros((len(y),2));wave=local.copy();report=dict(detected=False,period_MHz=1.09,selection='heldout_frequency_blocks',heldout_improvement_fraction=0.,amplitude_K=0.)
    if train.sum()<20 or test.sum()<20:return local,wave,report
    nu=(ff.mean(axis=1)-f.mean())*1000;x=nu/np.ptp(f*1000);nuisance=np.column_stack([np.ones(n),x,x*x]);mean=anchors.mean(axis=0)
    co=np.linalg.lstsq(nuisance[train],mean[train],rcond=None)[0];before=np.sqrt(np.mean((mean[test]-nuisance[test]@co)**2))
    trials=[]
    for period in np.linspace(1.04,1.14,41):
        # Exact finite-native-bin sinusoidal response.
        phase=2*np.pi*nu/period;df=float(np.median(np.diff(f)))*1000
        response=np.sinc(fit_bin*df/period)/np.sinc(df/period)
        cos=response*np.cos(phase);sin=response*np.sin(phase)
        design=np.column_stack([nuisance,cos,sin]);coeff=np.linalg.lstsq(design[train],mean[train],rcond=None)[0];prediction=design@coeff
        trials.append((float(np.mean((mean[train]-prediction[train])**2)),period,coeff,cos,sin,prediction))
    _,period,coeff,cos,sin,prediction=min(trials,key=lambda x:x[0]);after=np.sqrt(np.mean((mean[test]-prediction[test])**2));improvement=1-after/max(before,1e-15);amplitude=np.hypot(*coeff[-2:])
    detected=bool(improvement>=.1 and amplitude>=.03 and 1.04+1e-6<period<1.14-1e-6)
    report.update(detected=detected,period_MHz=float(period),heldout_improvement_fraction=float(improvement),amplitude_K=float(amplitude),heldout_RMS_before_K=float(before),heldout_RMS_after_K=float(after))
    if not detected:return local,wave,report
    design=np.stack([np.ones_like(v),v/limit,np.broadcast_to(cos,v.shape),np.broadcast_to(sin,v.shape)],axis=2);weights=good.astype(float)
    if np.any(good.sum(axis=1)<20):raise ValueError('Insufficient ripple fit anchors')
    for _ in range(8):
        matrix=np.einsum('rn,rni,rnj->rij',weights,design,design);target=np.einsum('rn,rni,rn->ri',weights,design,anchors)
        coeffs=np.linalg.solve(matrix,target[...,None])[...,0];res=anchors-np.einsum('rni,ri->rn',design,coeffs)
        sigma=np.array([max(.05,1.4826*np.median(abs(z[g]-np.median(z[g])))) for z,g in zip(res,good)])
        weights=good/np.sqrt(1+(res/sigma[:,None])**2)
    return coeffs[:,:2],coeffs[:,2:],report
