"""Explicit channel averaging, spectral smoothing, and backend gain transfer.

Frequencies may use any consistent unit. Data remain temperature or power
DENSITIES: normalized means, never sums. Channel edges assume a regular grid.
"""
from collections import OrderedDict
_OPERATOR_CACHE = OrderedDict()

import numpy as np
from scipy.sparse import csr_matrix
from scipy.ndimage import gaussian_filter1d


def channel_edges(centres):
    x = np.asarray(centres, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError('At least two finite channel centres are required')
    dx = np.diff(x)
    if np.any(dx <= 0) or not np.allclose(dx, dx[0], rtol=1e-6, atol=abs(dx[0])*1e-7):
        raise ValueError('Channel centres must be increasing and regularly spaced')
    return np.r_[x[0]-dx[0]/2, (x[:-1]+x[1:])/2, x[-1]+dx[-1]/2]


def overlap_rebin(data, input_freq, output_freq, axis=1, min_coverage=1.0):
    """Piecewise-constant channel-overlap averages; return values and coverage.

    Coverage is valid input bandwidth / requested output bandwidth. Incomplete
    output bins are NaN below min_coverage; no extrapolation or zero padding.
    Shifted neighbouring output bins can share input noise. This is a channel
    average model, not an exact instrumental-response or covariance correction.
    """
    if not 0 < min_coverage <= 1:
        raise ValueError('min_coverage must be in (0,1]')
    ie, oe = channel_edges(input_freq), channel_edges(output_freq)
    key = (len(input_freq), float(ie[0]), float(ie[-1]), len(output_freq), float(oe[0]), float(oe[-1]))
    operator = _OPERATOR_CACHE.get(key)
    if operator is None:
        rows, cols, weights = [], [], []
        for j, (lo, hi) in enumerate(zip(oe[:-1], oe[1:])):
            first = max(0, np.searchsorted(ie, lo, side='right')-1)
            last = min(len(ie)-1, np.searchsorted(ie, hi, side='left'))
            for i in range(first, last):
                w = max(0., min(hi, ie[i+1])-max(lo, ie[i]))/(hi-lo)
                if w:
                    rows.append(j); cols.append(i); weights.append(w)
        operator = csr_matrix((weights, (rows, cols)), shape=(len(output_freq), len(input_freq)))
        _OPERATOR_CACHE[key] = operator
        if len(_OPERATOR_CACHE) > 4:
            _OPERATOR_CACHE.popitem(last=False)
    a = np.moveaxis(np.asarray(data), axis, -1)
    if a.shape[-1] != len(input_freq):
        raise ValueError('Input frequency axis does not match data')
    flat = a.reshape(-1, a.shape[-1]); finite = np.isfinite(flat)
    support = (operator @ finite.astype(float).T).T
    summed = (operator @ np.where(finite, flat, 0.).T).T
    accepted = support >= min_coverage-1e-9
    result = np.divide(summed, support, out=np.full_like(summed, np.nan), where=accepted & (support > 0))
    shape = a.shape[:-1]+(len(output_freq),)
    return np.moveaxis(result.reshape(shape), -1, axis), np.moveaxis(support.reshape(shape), -1, axis)


def bounded_linear(data, input_freq, output_freq, axis=1):
    """Linear grid matching; NaN outside support or across an invalid pair."""
    x, q = np.asarray(input_freq), np.asarray(output_freq)
    if np.any(np.diff(x) <= 0):
        raise ValueError('Interpolation input must be increasing')
    a = np.moveaxis(np.asarray(data), axis, -1)
    i = np.clip(np.searchsorted(x, q, side='right')-1, 0, len(x)-2)
    w = (q-x[i])/(x[i+1]-x[i])
    lo, hi = a[..., i], a[..., i+1]
    valid = np.isfinite(lo) & np.isfinite(hi)
    out = np.where(valid, (1-w)*lo+w*hi, np.nan)
    # An exact match needs only its own channel, not a finite neighbour.
    out = np.where(w == 0, lo, out)
    out = np.where(w == 1, hi, out)
    out[..., (q < x[0]) | (q > x[-1])] = np.nan
    return np.moveaxis(out, -1, axis)


def gaussian_to_fwhm(data, n_g, native_fwhm_chan=1., axis=1, min_coverage=1.):
    """Smooth to a nominal Gaussian FWHM in native-channel units.

    Both input and target LSF are ASSUMED Gaussian. No change of sample grid.
    Subsequent interpolation contributes a further, non-Gaussian response.
    """
    if not np.isfinite(n_g) or not np.isfinite(native_fwhm_chan) or native_fwhm_chan <= 0 or n_g < native_fwhm_chan:
        raise ValueError('n_g must be >= positive native_fwhm_chan')
    if not 0 < min_coverage <= 1:
        raise ValueError('min_coverage must be in (0,1]')
    a = np.asarray(data)
    sigma = np.sqrt(n_g*n_g-native_fwhm_chan*native_fwhm_chan)/np.sqrt(8*np.log(2))
    if sigma < 1e-12:
        return a.copy()
    finite = np.isfinite(a)
    support = gaussian_filter1d(finite.astype(float), sigma, axis=axis, mode='constant', cval=0., truncate=4.)
    total = gaussian_filter1d(np.where(finite, a, 0.).astype(float), sigma, axis=axis, mode='constant', cval=0., truncate=4.)
    return np.divide(total, support, out=np.full_like(total, np.nan), where=support >= min_coverage-1e-9)


def frequency_gain(spec_data, spec_freq, tpsr, cal_freq, excluded_ranges=(), degree=3, clip_sigma=4., return_coefficients=False, evaluate=True):
    """Fit positive, smooth counts/K gains independently per time and auto-pol.

    Science power is averaged onto fully covered PSR channels. PSR channels
    touching excluded frequency intervals do not anchor the log-polynomial fit.
    This models a smooth relative bandpass, not narrow spectral structure.
    Additional polarizations receive NaN. Insufficient anchors fail explicitly.
    """
    if int(degree) != degree or degree < 0 or not np.isfinite(clip_sigma) or clip_sigma <= 0:
        raise ValueError('degree must be nonnegative integer and clip_sigma positive')
    spec_freq, cal_freq = np.asarray(spec_freq), np.asarray(cal_freq)
    all_edges = channel_edges(cal_freq)
    spec_edges = channel_edges(spec_freq)
    contained = (all_edges[:-1] >= spec_edges[0]) & (all_edges[1:] <= spec_edges[-1])
    if contained.sum() < degree+3:
        raise ValueError('Too little PSR/science band overlap')
    cal_freq = cal_freq[contained]
    tpsr = tpsr[:, :, contained]
    coarse, coverage = overlap_rebin(spec_data, spec_freq, cal_freq, axis=1)
    edges = channel_edges(cal_freq)
    allowed = np.ones(len(cal_freq), dtype=bool)
    for lo, hi in excluded_ranges:
        if hi < lo:
            raise ValueError('Excluded frequency intervals must increase')
        allowed &= ~((edges[:-1] <= hi) & (edges[1:] >= lo))
    midpoint = (spec_freq[0]+spec_freq[-1])/2; halfspan = (spec_freq[-1]-spec_freq[0])/2
    x = (cal_freq-midpoint)/halfspan; q = (spec_freq-midpoint)/halfspan
    result = np.full(spec_data.shape, np.nan, dtype=float) if evaluate else None; rejected = 0
    coefficients = np.full((len(spec_data),2,degree+1),np.nan) if return_coefficients else None
    for pol in range(min(2, spec_data.shape[2], tpsr.shape[1])):
        reference = tpsr[:, pol, :]
        ratio = np.divide(coarse[:, :, pol], reference, out=np.full_like(reference, np.nan, dtype=float), where=reference > 0)
        for row in range(len(spec_data)):
            good = allowed & np.isfinite(ratio[row]) & (ratio[row] > 0) & (coverage[row, :, pol] >= 1-1e-9)
            initial = int(good.sum())
            if initial < degree+3:
                raise ValueError(f'Insufficient positive calibration anchors at row {row}, pol {pol}')
            for iteration in range(4):
                coeff = np.polynomial.polynomial.polyfit(x[good], np.log(ratio[row, good]), degree)
                residual = np.log(ratio[row, good])-np.polynomial.polynomial.polyval(x[good], coeff)
                centre = np.median(residual); scatter = 1.4826*np.median(np.abs(residual-centre))
                if iteration == 3 or scatter < 1e-12:
                    break
                indices = np.flatnonzero(good); keep = np.abs(residual-centre) <= clip_sigma*scatter
                if keep.all() or keep.sum() < degree+3:
                    break
                good[indices[~keep]] = False
            # Require interpolation support on BOTH sides of the HI mask.
            if np.min(x[good]) >= 0 or np.max(x[good]) <= 0:
                raise ValueError('Calibration anchors must span the science band centre')
            if return_coefficients:
                coefficients[row,pol]=coeff
            if evaluate:
                gain = np.exp(np.polynomial.polynomial.polyval(q, coeff))
                # No extrapolation beyond the calibration reference channel support.
                supported = (spec_freq >= all_edges[0]) & (spec_freq <= all_edges[-1])
                result[row, supported, pol] = gain[supported]
            rejected += initial-int(good.sum())
    info = {'degree':int(degree), 'sigma_clip':float(clip_sigma), 'rejected_anchors':rejected,
            'model':'log-polynomial relative backend gain; top-hat PSR matching assumption'}
    if return_coefficients:
        info['coefficients']=coefficients
        info['reference_edges_ghz']=[float(all_edges[0]),float(all_edges[-1])]
    return result, info


def _affine_arrays(data, frequency_scale, axis):
    """Flatten non-frequency dimensions, broadcasting time-axis factors."""
    array=np.asarray(data);axis=axis % array.ndim
    if array.ndim < 2 or axis == 0:
        raise ValueError('Time axis must be first; frequency axis must be different')
    moved=np.moveaxis(array,axis,-1)
    scale=np.asarray(frequency_scale,dtype=float)
    if scale.ndim==0:scale=np.full(array.shape[0],scale)
    if scale.shape != (array.shape[0],) or not np.isfinite(scale).all() or np.any(scale<=0):
        raise ValueError('One positive finite frequency factor per spectrum required')
    leading=moved.shape[:-1]
    scales=np.broadcast_to(scale.reshape((len(scale),)+(1,)*(len(leading)-1)),leading).reshape(-1)
    return moved.reshape(-1,moved.shape[-1]),scales,moved.shape,axis


def affine_overlap_rebin(data,input_freq,output_freq,frequency_scale,axis=1,min_coverage=1.):
    """Per-spectrum Doppler alignment and bin means in one vectorized pass.

    Cumulative integrals of native piecewise-constant channels avoid a Python
    loop or a sparse matrix per spectrum. Factors multiply native frequencies.
    Missing support is propagated through a separate cumulative valid-width sum.
    """
    if not 0 < min_coverage <= 1:raise ValueError('min_coverage must be in (0,1]')
    ie,oe=channel_edges(input_freq),channel_edges(output_freq)
    flat,scales,shape,axis=_affine_arrays(data,frequency_scale,axis)
    if flat.shape[-1]!=len(input_freq):raise ValueError('Frequency axis length mismatch')
    dx=float((input_freq[-1]-input_freq[0])/(len(input_freq)-1));n=len(input_freq)
    positions=(oe[None,:]/scales[:,None]-ie[0])/dx
    clipped=np.clip(positions,0,n);indices=np.floor(clipped).astype(np.int64)
    fraction=clipped-indices;sample=np.minimum(indices,n-1)
    valid=np.isfinite(flat);values=np.where(valid,flat,0.)
    sums=np.pad(np.cumsum(values,axis=1,dtype=float),((0,0),(1,0)))
    counts=np.pad(np.cumsum(valid,axis=1,dtype=float),((0,0),(1,0)))
    integral=np.take_along_axis(sums,indices,axis=1)+fraction*np.take_along_axis(values,sample,axis=1)
    support=np.take_along_axis(counts,indices,axis=1)+fraction*np.take_along_axis(valid,sample,axis=1)
    integrated=np.diff(integral,axis=1);width=np.diff(support,axis=1)
    coverage=width/np.diff(positions,axis=1)
    result=np.divide(integrated,width,out=np.full_like(integrated,np.nan),where=(width>0)&(coverage>=min_coverage-1e-9))
    output_shape=shape[:-1]+(len(output_freq),)
    return np.moveaxis(result.reshape(output_shape),-1,axis),np.moveaxis(coverage.reshape(output_shape),-1,axis)


def affine_linear(data,input_freq,output_freq,frequency_scale,axis=1):
    """Per-spectrum bounded linear grid matching after Gaussian smoothing."""
    channel_edges(input_freq)
    flat,scales,shape,axis=_affine_arrays(data,frequency_scale,axis)
    q=np.asarray(output_freq)[None,:]/scales[:,None]
    positions=(q-input_freq[0])/((input_freq[-1]-input_freq[0])/(len(input_freq)-1))
    indices=np.clip(np.floor(positions).astype(np.int64),0,len(input_freq)-2)
    fraction=positions-indices
    lo=np.take_along_axis(flat,indices,axis=1);hi=np.take_along_axis(flat,indices+1,axis=1)
    out=np.where(np.isfinite(lo)&np.isfinite(hi),(1-fraction)*lo+fraction*hi,np.nan)
    out=np.where(fraction==0,lo,out);out=np.where(fraction==1,hi,out)
    out[(q<input_freq[0])|(q>input_freq[-1])]=np.nan
    return np.moveaxis(out.reshape(shape[:-1]+(len(output_freq),)),-1,axis)
