"""Physical-axis diagnostics for FISH; display binning never changes science arrays."""
from pathlib import Path
import warnings
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

COLORS = ['#1674a5', '#d97720']
POL_NAMES = ['AA', 'BB']
STYLE = {'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
         'axes.spines.top': False, 'axes.spines.right': False,
         'savefig.dpi': 160, 'figure.facecolor': 'white'}


def _median(values, axis):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmedian(np.where(np.isfinite(values), values, np.nan), axis=axis)


def shared_limits(arrays, percentiles=(2, 98)):
    samples = []
    for array in arrays:
        a = np.asarray(array).ravel()
        a = a[::max(1, len(a)//200000)]
        samples.append(a[np.isfinite(a)])
    finite = np.concatenate(samples)
    if not len(finite):
        return (0, 1)
    lo, hi = np.percentile(finite, percentiles)
    if hi <= lo:
        hi = lo + max(abs(lo)*0.01, 1e-6)
    return float(lo), float(hi)


def _edges(centers):
    centers = np.asarray(centers)
    if len(centers) == 1:
        return np.array([centers[0]-.5, centers[0]+.5])
    middle = (centers[:-1]+centers[1:])/2
    return np.r_[centers[0]-(middle[0]-centers[0]), middle, centers[-1]+centers[-1]-middle[-1]]


def _display_bins(data, time, freq, max_time=512, max_freq=1200):
    """Block-mean only the displayed image; include partial final blocks."""
    bt = max(1, int(np.ceil(len(time)/max_time)))
    bf = max(1, int(np.ceil(len(freq)/max_freq)))
    finite = np.isfinite(data)
    values = np.where(finite, data, 0.).astype(float)
    counts = finite.astype(np.int32)
    for size, axis in [(bf, 1), (bt, 0)]:
        indices = np.arange(0, values.shape[axis], size)
        values = np.add.reduceat(values, indices, axis=axis)
        counts = np.add.reduceat(counts, indices, axis=axis)
    image = np.divide(values, counts, out=np.full_like(values, np.nan), where=counts>0)
    t = np.array([np.mean(time[i:i+bt]) for i in range(0,len(time),bt)])
    f = np.array([np.mean(freq[i:i+bf]) for i in range(0,len(freq),bf)])
    return image, t, f, bt, bf


def _save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(stem)+'.pdf')
    fig.savefig(str(stem)+'.png', dpi=160)
    plt.close(fig)


def waterfall(data, mjd, freq_mhz, stem, title, pol=0, unit='K',
              reference=None, reference_label='Smoothed', limits=None, note=''):
    data = np.asarray(data)
    mjd, freq_mhz = np.asarray(mjd), np.asarray(freq_mhz)
    if data.shape != (len(mjd), len(freq_mhz)):
        raise ValueError('Diagnostic axes do not match data shape')
    time = (mjd-mjd[0])*1440
    order = np.argsort(freq_mhz)
    freq_mhz, data = freq_mhz[order], data[:, order]
    if reference is not None:
        reference = np.asarray(reference)[:, order]
    image, t, f, bt, bf = _display_bins(data, time, freq_mhz)
    lo, hi = limits or shared_limits([data])
    dt = np.median(np.diff(mjd))*86400 if len(mjd)>1 else np.nan
    df = np.median(np.diff(freq_mhz))*1000 if len(freq_mhz)>1 else np.nan
    finite = np.isfinite(data)
    clipped = np.count_nonzero(finite & ((data<lo)|(data>hi)))/max(1,np.count_nonzero(finite))*100
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(11.8, 7.6), constrained_layout=True)
        gs = fig.add_gridspec(2, 3, width_ratios=[4.6, 1.45, .14], height_ratios=[3.3, 1.2])
        ax = fig.add_subplot(gs[0,0]); side = fig.add_subplot(gs[0,1], sharey=ax)
        spectrum = fig.add_subplot(gs[1,0], sharex=ax); info = fig.add_subplot(gs[1,1:]);info.axis('off')
        cmap = plt.get_cmap('viridis').copy();cmap.set_bad('#d9dde2')
        mesh = ax.pcolormesh(_edges(f), _edges(t), image, cmap=cmap, vmin=lo, vmax=hi,
                             shading='flat', rasterized=True)
        fig.colorbar(mesh,cax=fig.add_subplot(gs[0,2]),label=unit,extend='both')
        spectrum.plot(freq_mhz, _median(data,0), color=COLORS[pol], lw=1, label='Native-data median')
        side.plot(_median(data,1), time, color=COLORS[pol], lw=1, label='Native-data median')
        if reference is not None:
            spectrum.plot(freq_mhz, _median(reference,0), color='#bd422d', lw=1.4, label=reference_label)
            side.plot(_median(reference,1), time, color='#bd422d', lw=1.4)
        # Make missing-time intervals visible rather than implying continuous coverage.
        if len(time)>1:
            gap = np.where(np.diff(time)>1.5*np.median(np.diff(time)))[0]
            for i in gap:
                ax.axhspan(time[i]+dt/120, time[i+1]-dt/120, color='#d9dde2', zorder=3)
        ax.tick_params(labelbottom=False);side.tick_params(labelleft=False)
        ax.set_ylabel('Elapsed time (min)');side.set_xlabel('Band median\n'+unit)
        spectrum.set(xlabel='Frequency (MHz)',ylabel='Time median\n'+unit)
        spectrum.legend(loc='best',fontsize=8,frameon=False)
        spectrum.grid(alpha=.18);side.grid(alpha=.18)
        ax.set_xlim(freq_mhz[0],freq_mhz[-1]);ax.set_ylim(_edges(t)[0],_edges(t)[-1])
        fig.suptitle(f'{title} — {POL_NAMES[pol]}\nUTC start: {Time(mjd[0],format="mjd").isot}',fontsize=13)
        statistics = (f'{data.shape[0]:,} × {data.shape[1]:,} native samples\n'
                      f'Δt: {dt:.6f} s | Δν: {df:.3f} kHz\n'
                      f'Finite: {100*finite.mean():.3f}% | color-clipped: {clipped:.2f}%\n'
                      f'Color limits: {lo:.4g}–{hi:.4g} (shared AA/BB)\n'
                      f'Display bins: {bt} time × {bf} frequency\n'
                      'Marginals: full-resolution arrays\n'
                      'Clipped colors do not remove data.\n')
        if note:
            statistics += textwrap.fill(note, width=45)
        info.text(0,.98,statistics,va='top',fontsize=8.2,linespacing=1.35)
        _save(fig,stem)


def tcal_reference(freq_mhz, native, stem, target_freq_mhz=None, processed=None, title='Noise-diode temperature reference'):
    with plt.rc_context(STYLE):
        rows=3 if processed is not None else 2
        fig, axes=plt.subplots(rows,1,figsize=(10.8,7.5 if rows==3 else 5.7),sharex=True,constrained_layout=True)
        for p in range(2):
            axes[p].plot(freq_mhz,native[p],color='#9ba4ad',lw=.7,label='Native reference')
            if processed is not None:
                axes[p].plot(target_freq_mhz,processed[p],color=COLORS[p],lw=1.15,label='Smoothed and interpolated')
                unfiltered=np.interp(target_freq_mhz,freq_mhz,native[p])
                supported=(target_freq_mhz>=np.min(freq_mhz))&(target_freq_mhz<=np.max(freq_mhz))
                axes[2].plot(target_freq_mhz[supported],1000*(processed[p,supported]-unfiltered[supported]),color=COLORS[p],lw=.9,label=POL_NAMES[p])
            axes[p].set_ylabel(f'{POL_NAMES[p]} Tcal (K)');axes[p].grid(alpha=.18);axes[p].legend(frameon=False,fontsize=8)
        if rows==3:
            axes[2].axhline(0,color='#65717d',lw=.7);axes[2].set_ylabel('Processing difference\n(mK)');axes[2].legend(frameon=False,fontsize=8)
            axes[2].grid(alpha=.18)
        axes[-1].set_xlabel('Frequency (MHz)');fig.suptitle(title)
        _save(fig,stem)


def scale_transfer(mjd, spec_median, tpsr_median, factor, stem, title):
    time=(np.asarray(mjd)-mjd[0])*1440
    factor=np.asarray(factor)
    if factor.ndim==3:factor=factor[:,0,:]
    with plt.rc_context(STYLE):
        fig,axes=plt.subplots(3,1,figsize=(10.8,7.6),sharex=True,constrained_layout=True)
        for p in range(2):
            axes[0].plot(time,tpsr_median[:,p],color=COLORS[p],lw=1,label=POL_NAMES[p])
            axes[1].plot(time,spec_median[:,p],color=COLORS[p],lw=1)
            normalization=_median(factor[:,p],None)
            axes[2].plot(time,100*(factor[:,p]/normalization-1),color=COLORS[p],lw=1)
        axes[0].set_ylabel('PSR band median (K)');axes[1].set_ylabel('Spectral band median\n(arb. units)')
        axes[2].set_ylabel('Conversion factor\nrelative to median (%)');axes[2].set_xlabel('Elapsed time (min)')
        axes[0].legend(frameon=False,ncol=2);axes[2].axhline(0,color='#65717d',lw=.7)
        for ax in axes:ax.grid(alpha=.18)
        fig.suptitle(f'{title}\nUTC start: {Time(mjd[0],format="mjd").isot}',fontsize=12)
        _save(fig,stem)
