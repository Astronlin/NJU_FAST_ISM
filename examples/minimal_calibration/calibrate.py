#!/usr/bin/env python3
"""Real Sun17_A_01/M01 calibration example using production FISH functions."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from fastfish.calibration.extracted import rebuild_reference, extracted_band
from fastfish.calibration.reference import smoothing_rows
from fastfish.processing.preview import calibrate

HERE = Path(__file__).resolve().parent
REST_GHZ = 1.420405751768


def load_inputs():
    data = HERE / 'data'
    provenance = json.loads((data / 'provenance.json').read_text())
    for name, record in provenance['files'].items():
        content = (data / name).read_bytes()
        if len(content) != record['bytes'] or hashlib.sha256(content).hexdigest() != record['sha256']:
            raise ValueError('Input checksum mismatch: ' + name)
    with fits.open(data / 'science.fits', memmap=True) as h:
        tab = h[1].data
        raw = tab['DATA'].astype(float)
        frequency = (float(tab['FREQ'][0]) + np.arange(raw.shape[1]) * float(tab['CHAN_BW'][0])) / 1000
        # Original SCI UTOBS is integration start; CAL times are already midpoints.
        times = tab['UTOBS'].copy() + tab['EXPOSURE'].copy() / (2 * 86400)
    return raw, frequency, times, provenance


def build_reference(frequency, times, n_t_cal=10, n_t_apply=10):
    for value in (n_t_cal, n_t_apply):
        smoothing_rows(value)
        if value > 10:
            raise ValueError('This compact extract supports smoothing windows up to 10; larger windows need more CAL context')
    return rebuild_reference(
        [HERE / 'data/on_off.fits'], float(times[0]), float(times[-1]), frequency,
        HERE / 'data/tcal_reference.fits', 'M01', n_t=n_t_apply, n_t_cal=n_t_cal)


def apply_reference(raw, frequency, times, pack):
    tpsr = interp1d(pack['mjd'], pack['smoothed'], axis=0, bounds_error=True)(times)
    cfg = dict(rest_freq_ghz=REST_GHZ, gain_degree=3, protect_kms=250., rfi_catalog='packaged')
    # Deliberately no pointing/Doppler calculation: scale=1 means topocentric.
    # This runs the production chunk RFI search and frequency-gain estimator.
    result = calibrate(raw, frequency, tpsr, pack['frequency_ghz'],
                       np.ones(len(times)), [], cfg)
    return result, tpsr


def run(output, n_t_cal=10, n_t_apply=10):
    output = Path(output).resolve()
    if output == HERE or HERE in output.parents:
        raise ValueError('Choose an output directory outside the released example')
    raw, frequency, times, provenance = load_inputs()
    pack = build_reference(frequency, times, n_t_cal, n_t_apply)
    result, tpsr = apply_reference(raw, frequency, times, pack)
    temperature = result['pol']
    expected = json.loads((HERE / 'expected/regression_reference.json').read_text())
    indices = expected['channel_indices']
    reference = np.asarray(expected['temperature_K'])
    error = temperature[:, indices, :] - reference
    default = n_t_cal == 10 and n_t_apply == 10
    # This is repeatability against a frozen production calculation, NOT truth.
    if default:
        np.testing.assert_allclose(temperature[:, indices, :], reference,
                                   atol=expected['tolerance_atol_K'], rtol=expected['tolerance_rtol'])
        np.testing.assert_allclose(result['coefficients'], expected['coefficients'], atol=1e-8, rtol=1e-8)
    assert np.isfinite(temperature).all(), 'Unexpected missing calibrated samples'
    # Real CAL can have nonpositive values in excluded H I channels. The
    # production gain fit requires positive supported anchors, not every bin.
    summary = dict(dataset='Sun17_A_01 real observation', beam='M01',
        science_shape=list(raw.shape), reference_shape=list(pack['raw'].shape),
        n_t_cal=n_t_cal, n_t_apply=n_t_apply,
        science_duration_s=provenance['science_duration_s'],
        midpoint_MJD_range=times[[0,-1]].tolist(),
        finite_fraction=float(np.isfinite(temperature).mean()),
        median_Tpsr_K=np.median(tpsr,axis=(0,2)).tolist(),
        nonpositive_reference_fraction_by_pol=np.mean(tpsr<=0,axis=(0,2)).tolist(),
        nonfinite_reference_fraction_by_pol=np.mean(~np.isfinite(tpsr),axis=(0,2)).tolist(),
        nonfinite_reference_frequencies_MHz=(pack['frequency_ghz'][np.any(~np.isfinite(tpsr),axis=(0,1))]*1000).tolist(),
        nonpositive_reference_frequencies_MHz=(pack['frequency_ghz'][np.any(tpsr<=0,axis=(0,1))]*1000).tolist(),
        median_calibrated_K=np.median(temperature,axis=(0,1)).tolist(),
        detected_RFI_intervals_MHz=result['intervals'],
        regression_compared=default,
        maximum_sample_change_from_default_K=float(np.max(abs(error))),
        cropped_vs_original_max_K=provenance['calibrated_comparison_max_K'],
        gain_coefficients_shape=list(result['coefficients'].shape),
        checks_passed=True,
        note='Real-data numerical repeatability; no known astrophysical truth, no absolute calibration-accuracy claim. No baseline subtraction or Doppler correction.')
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'calibrated.npz', frequency_ghz=frequency,
                        midpoint_mjd=times, temperature_K=temperature,
                        average_K=result['average'], gain_coefficients=result['coefficients'])
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    plot_diagnostics(raw, frequency, times, pack, result, output / 'calibration.png')
    print(json.dumps(summary, indent=2))
    print('Outputs:', output)
    return summary


def plot_diagnostics(raw, f, times, pack, result, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tt, on, off, cf = extracted_band(str(HERE / 'data/on_off.fits'), float(f[0]), float(f[-1]))
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    velocity = (1-f/REST_GHZ)*299792.458
    # Display-only means retain native spectra for all calibration operations.
    display_bin = 32
    def display_mean(values):
        return np.asarray(values).reshape(-1,display_bin).mean(axis=1)
    fm = display_mean(f)*1000; vm = display_mean(velocity)
    for p, label in enumerate(('AA','BB')):
        color=('#1764ab','#c24f17')[p]
        ax[0,0].plot((tt-tt[0])*86400,np.median((on-off)/20480*(20480+28672)/2,axis=2)[:,p],color=color,alpha=.35)
        from scipy.ndimage import uniform_filter1d
        sm=uniform_filter1d((on-off)/20480*(20480+28672)/2,pack['info']['onoff_smoothing_rows'],axis=0,mode='nearest')
        ax[0,0].plot((tt-tt[0])*86400,np.median(sm,axis=2)[:,p],color=color,label=label)
        ax[0,1].plot((pack['mjd']-tt[0])*86400,np.median(pack['raw'],axis=2)[:,p],color=color,alpha=.35)
        ax[0,1].plot((pack['mjd']-tt[0])*86400,np.median(pack['smoothed'],axis=2)[:,p],color=color,label=label)
        gain=np.divide(raw[0,:,p],result['pol'][0,:,p])
        ax[1,0].plot(fm,display_mean(gain),color=color,label=label)
        ax[1,1].plot(vm,display_mean(result['pol'][0,:,p]),color=color,label=label)
    ax[0,0].set(xlabel='Time since first CAL midpoint (s)',ylabel='Duty-adjusted ON−OFF (counts)',title='Sun17_A_01 M01 | Diode increment')
    ax[0,1].set(xlabel='Time since first CAL midpoint (s)',ylabel='Band-median Tpsr (K)',title='Temperature reference')
    ax[1,0].set(xlabel='Native frequency (MHz)',ylabel='Gain (counts/K)',title='Frequency-dependent gain | first SCI row')
    ax[1,1].set(xlabel='Topocentric radio velocity (km/s)',ylabel='Temperature (K)',xlim=(-150,150),title='Calibrated AA/BB | no baseline subtraction')
    for a in ax.flat:a.legend();a.grid(alpha=.2)
    fig.savefig(path,dpi=150);plt.close(fig)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--n-t-cal',type=int,default=10)
    parser.add_argument('--n-t-apply',type=int,default=10)
    args=parser.parse_args()
    run(args.output,args.n_t_cal,args.n_t_apply)
