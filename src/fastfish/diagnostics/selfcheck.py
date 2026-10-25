"""Synthetic gain, baseline, protected-line and rebin regression."""
from pathlib import Path
import json
import numpy as np
from fastfish.processing.spectral import frequency_gain,affine_overlap_rebin
from fastfish.rfi.detection import fit_linear_anchors,frequency_mask,C
from fastfish.baseline.broad import fit_global_broad
f=1.4047746658325195+np.arange(65536)*476.837158203125/1e9
rest=1.420405751768;scale=np.array([1.0001043]);v=(1-f*scale[0]/rest)*C
line=100*np.exp(-.5*((v-7)/10)**2);temperature=40+.1*v/1000+line
intervals=[(1405.3,1406.3),(1417.3,1418.4),(1421.4,1422.6),(1433.3,1434.5)]
# RFI is compactly supported inside the known masks in this controlled test.
rfi=np.zeros(len(f))
for lo,hi in intervals:rfi[(f*1000>lo+.1)&(f*1000<hi-.1)]=2.
x=(f-f.mean())/np.ptp(f);gain=np.exp(np.stack([np.log(500)+.15*x+.03*x*x,np.log(400)-.1*x+.02*x*x],axis=1))
raw=(temperature+rfi)[None,:,None]*gain[None]
cf=np.arange(f[0]-.0002,f[-1]+.0002,.000125)
# Reference contains the true smooth temperature continuum on its own grid.
cv=(1-cf*scale[0]/rest)*C;tpsr=np.broadcast_to((40+.1*cv/1000)[None,None,:],(1,2,len(cf)))
exclude=[((1-700/C)*rest,(1+700/C)*rest)]+[(a/1000,b/1000) for a,b in intervals]
estimated,_=frequency_gain(raw,f,tpsr,cf,excluded_ranges=exclude,degree=3)
cal=np.mean(raw/estimated,axis=2)
global_model,gco,global_report=fit_global_broad(cal,f,scale,rest,intervals)
co,info=fit_linear_anchors(cal-global_model,f,scale,rest,intervals,fit_limit=1500.)
residual=cal-global_model-co[:,0,None]-co[:,1,None]*v[None,:]/1500
# Identified RFI remains in the native science spectrum; masks affect fits only.
assert np.isfinite(residual).all()
assert np.median(residual[0,(f*1000>1405.5)&(f*1000<1406.1)])>1.9
q=np.arange((1-100/C)*rest,(1+100/C)*rest,2*476.837158203125/1e9)
out,cov=affine_overlap_rebin(residual,f,q,scale,axis=1)
native_dv=476.837158203125/1e9*scale[0]/rest*C;output_dv=2*476.837158203125/1e9/rest*C
expected=line.sum()*native_dv;area=np.nansum(out)*output_dv;relative=abs(area/expected-1)
# This challenge includes a full-band slope, not an exact sinusoidal baseline.
# Require <0.01% line-area bias; record the actual model mismatch, not zero.
assert relative<1e-4,(area,expected,relative)
# An explicitly invalid native channel must make intersecting output bins invalid.
j=np.argmin(abs(v));bad=residual.copy();bad[0,j]=np.nan
changed,support=affine_overlap_rebin(bad,f,q,scale,axis=1)
assert np.sum(~np.isfinite(changed))>=1
assert np.all(support[~np.isfinite(changed)]<1-1e-9)
report={'known_line_area_K_kms':float(expected),'recovered_area_K_kms':float(area),'relative_area_error':float(relative),
    'invalid_native_channel_propagated':True,'identified_RFI_retained_in_science':True,'fit_anchors_each_side':[info['minimum_left_anchors'],info['minimum_right_anchors']],
    'scope':'Noise-free gain, broad-plus-local baseline and line-area challenge with known masks. The true continuum includes a full-band slope, so the finite-period broad model is imperfect; tolerance is 0.01 percent line-area bias. Does not validate real RFI classifications, absolute calibration, or baseline truth.'}
# Exact requested model: known broad sinusoid plus a protected narrow line.
phase=2*np.pi*(f-f.mean())*1000/23.
true_broad=40+.4*np.cos(phase)+.2*np.sin(phase)
signal=true_broad[None]+line[None]
model,gcoeff,details=fit_global_broad(signal,f,scale,rest,[])
lcoeff,_=fit_linear_anchors(signal-model,f,scale,rest,[],fit_limit=1500.)
remaining=signal[0]-model[0]-lcoeff[0,0]-lcoeff[0,1]*v/1500
max_error=float(np.max(abs(remaining-line)))
assert max_error<1e-8,max_error
changed=signal.copy();protected=abs(v)<250;changed[0,protected]+=1000*np.sin(v[protected])
model_changed,_,_=fit_global_broad(changed,f,scale,rest,[])
assert np.array_equal(model,model_changed)
report['exact_broad_model_test']={'known_period_MHz':23.,'recovered_period_MHz':details['period_MHz'],'maximum_line_error_K':max_error,'protected_fit_invariance_with_fixed_masks':True}
print(json.dumps(report,indent=2))
