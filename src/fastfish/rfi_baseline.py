"""Conservative harmonic RFI tracking and protected native linear baselines.

Detection averages a short time chunk and 128 native frequency channels ONLY
for mask estimation. Fits average polynomial bases exactly; science spectra
are never smoothed here. No feature is subtracted as an RFI model.
"""
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.signal import find_peaks
from .rfi_catalog import load_catalog,load_published_statistics,catalog_path,catalog_fingerprint
from functools import lru_cache
C=299792.458

def highpass(y):
    return gaussian_filter1d(y-median_filter(y,49,mode='nearest'),1.,mode='nearest')

def learn_priors(temperature,frequency_ghz,scales,labels,rest,pattern=None):
    if pattern is None:
        pattern=load_catalog()['patterns']['repeated_broad_features']
    spacings=pattern['spacing_MHz'];tolerance=pattern['match_tolerance_MHz']
    f=np.asarray(frequency_ghz)*1000;df=float(np.median(np.diff(f)));result={}
    for scan in sorted(set(s.split('/')[0] for s in labels)):
        ids=np.array([s.startswith(scan+'/') for s in labels])
        yy=np.median(temperature[ids].mean(axis=2),axis=0);z=highpass(yy)
        v=(1-f/1000*np.median(scales[ids])/rest)*C
        peaks,_=find_peaks(z,prominence=.14,distance=4)
        peaks=[j for j in peaks if z[j]>.15 and abs(v[j])>120]
        selected=[]
        for j in peaks:
            partners=[k for k in peaks if k!=j and min(abs(abs(f[k]-f[j])-spacing) for spacing in spacings)<tolerance]
            if partners:selected.append({'centre_MHz':float(f[j]),'median_excess_K':float(z[j]),
                'partners_MHz':[float(f[k]) for k in partners]})
        result[scan]=selected
    return result

@lru_cache(maxsize=4)
def _published_candidates(path, fingerprint):
    catalog=load_catalog(path)
    if 'published_statistics' not in catalog:return ()
    table=load_published_statistics(catalog,Path(path).parent)
    return tuple(row[0] for row in table['rows'] if row[2]>0 or row[4]>0)


def track_intervals(raw_anchor_mean,freq_anchor_ghz,priors,scale,rest,temperature_scale=40.,protect=250.,catalog=None):
    """Prior-guided detection, not prior-driven masking.

    Require a positive observed highpass peak above both 0.15 K-equivalent
    and six robust scatter estimates. Published 0.1-MHz bins provide search
    locations only. Nondetections and protected H I never generate intervals.
    Chunk-averaged detection can miss short-lived, weak or uncatalogued RFI.
    """
    f=np.asarray(freq_anchor_ghz)*1000;y=np.asarray(raw_anchor_mean,dtype=float)
    if not np.isfinite(y).all() or np.median(y)<=0:raise ValueError('Invalid RFI detection spectrum')
    normalized=y/np.median(y)*temperature_scale;z=highpass(normalized)
    velocity=(1-f/1000*scale/rest)*C;core=np.abs(velocity)<=protect
    clean=z[~core];scatter=1.4826*np.median(abs(clean-np.median(clean))) if len(clean) else 0.
    detection_limit=max(.15,6*float(scatter))
    candidates=[dict(p,search_halfwidth_MHz=.8,source='observation_profile') for p in priors]
    if catalog is not None:
        path=catalog_path(catalog)
        candidates.extend(dict(centre_MHz=x,search_halfwidth_MHz=.1,source='published_statistics')
            for x in _published_candidates(str(path),catalog_fingerprint(path)) if f[0]-.1<=x<=f[-1]+.1)
    intervals=[];records=[];seen=set()
    for prior in candidates:
        centre=prior['centre_MHz'];width=prior['search_halfwidth_MHz']
        search=np.flatnonzero((abs(f-centre)<=width)&~core)
        if not len(search):continue
        j=search[np.argmax(z[search])];peak=float(z[j])
        if width>.1 and abs(f[j]-centre)>.7:
            j=int(np.argmin(abs(f-centre)));peak=float(z[j])
        detected=bool(not core[j] and peak>=detection_limit)
        record=dict(prior_MHz=centre,tracked_MHz=float(f[j]),raw_equivalent_peak_K=peak,
            detection_threshold_K_equivalent=detection_limit,detected=detected,source=prior['source'],interval_MHz=None,weak_prior_mask=False)
        records.append(record)
        if not detected or j in seen:continue
        seen.add(j);left=right=j;threshold=max(.05,.025*peak)
        while left>0 and not core[left-1] and z[left-1]>threshold and f[j]-f[left-1]<1.25:left-=1
        while right<len(f)-1 and not core[right+1] and z[right+1]>threshold and f[right+1]-f[j]<1.25:right+=1
        lo=min(float(f[left]-.15),float(f[j]-.35));hi=max(float(f[right]+.15),float(f[j]+.35))
        core_freq=f[core]
        if len(core_freq):
            half=float(np.median(np.diff(f)))/2
            if f[j]<core_freq[0]:hi=min(hi,float(core_freq[0]-half))
            else:lo=max(lo,float(core_freq[-1]+half))
        intervals.append((lo,hi));record['interval_MHz']=[lo,hi]
    intervals.sort();merged=[]
    for lo,hi in intervals:
        if merged and lo<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],hi)
        else:merged.append([lo,hi])
    return merged,records

def frequency_mask(freq_ghz,intervals,half_channel=0.):
    f=np.asarray(freq_ghz)*1000;mask=np.zeros(f.shape,dtype=bool)
    for lo,hi in intervals:mask|=(f+half_channel>=lo)&(f-half_channel<=hi)
    return mask

def fit_linear_anchors(data,freq_ghz,scales,rest,intervals,fit_bin=128,fit_limit=1000.,protect=250.,extra_exclusions=(),allow_missing=False):
    """Fit native-channel means, excluding RFI and ALL protected input values.

    Return c0,c1 for x=v_LSRK/fit_limit. Soft-L1 IRLS limits remaining outliers.
    At least eight whole anchors on each side are required; fail explicitly.
    """
    a=np.asarray(data,dtype=float);f=np.asarray(freq_ghz);scales=np.asarray(scales)
    n=len(f)//fit_bin;stop=n*fit_bin;ff=f[:stop].reshape(n,fit_bin)
    anchor=a[:,:stop].reshape(len(a),n,fit_bin).mean(axis=2)
    velocity=(1-ff.mean(axis=1)[None,:]*scales[:,None]/rest)*C
    vl=(1-ff[:,0][None,:]*scales[:,None]/rest)*C
    vr=(1-ff[:,-1][None,:]*scales[:,None]/rest)*C
    rfi=frequency_mask(f,intervals).reshape(n,fit_bin).any(axis=1)
    good=(vl<=fit_limit)&(vr>=-fit_limit)&((vr>protect)|(vl<-protect))&~rfi[None,:]&np.isfinite(anchor)
    for lo,hi in extra_exclusions:good&=~((vl>=lo)&(vr<=hi))
    left=(good&(velocity<-protect)).sum(axis=1);right=(good&(velocity>protect)).sum(axis=1)
    supported=(left>=8)&(right>=8)
    if not np.all(supported):
        if not allow_missing:
            raise ValueError('Insufficient clean baseline support: '+str([left.min(),right.min()]))
        # Diagnostic-only opt-in: never weaken the science fit's support rule.
        coeff=np.full((len(a),2),np.nan)
        if np.any(supported):
            coeff[supported],_=fit_linear_anchors(a[supported],f,scales[supported],rest,intervals,
                fit_bin,fit_limit,protect,extra_exclusions)
        return coeff,dict(minimum_left_anchors=int(left.min()),minimum_right_anchors=int(right.min()),
            diagnostic_available=supported,left_anchors=left,right_anchors=right)

    x=velocity/fit_limit;w=good.astype(float);coeff=np.empty((len(a),2))
    for iteration in range(8):
        s0=w.sum(axis=1);s1=(w*x).sum(axis=1);s2=(w*x*x).sum(axis=1)
        t0=(w*anchor).sum(axis=1);t1=(w*x*anchor).sum(axis=1);det=s0*s2-s1*s1
        if np.any(det<=1e-10):raise ValueError('Singular masked linear fit')
        coeff[:,0]=(t0*s2-t1*s1)/det;coeff[:,1]=(s0*t1-s1*t0)/det
        residual=anchor-coeff[:,0,None]-coeff[:,1,None]*x
        scatter=np.array([max(.05,1.4826*np.median(np.abs(rr[gg]-np.median(rr[gg])))) for rr,gg in zip(residual,good)])
        w=good/np.sqrt(1+(residual/scatter[:,None])**2)
    info={'minimum_left_anchors':int(left.min()),'minimum_right_anchors':int(right.min()),
        'anchor_RMS_K':np.array([np.sqrt(np.mean(rr[gg]**2)) for rr,gg in zip(residual,good)]),
        'anchor_used':good,'anchor_frequency_ghz':ff.mean(axis=1),'anchor_residual_K':residual}
    return coeff,info
