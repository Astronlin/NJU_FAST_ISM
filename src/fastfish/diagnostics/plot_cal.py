"""Standalone versions of the same FISH calibration diagnostics."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from fastfish.diagnostics.calibration import waterfall, shared_limits, tcal_reference


def plot_Tcal(Tcalfits='Tcal/CAL.20201014.low.W.fits', beam=1, outdir='.'):
    with fits.open(Tcalfits) as hdus:
        frequency = hdus[1].data['FREQ'][0].copy()
        temperatures = hdus[1].data['TCAL'][0, beam-1].T.copy()
    tcal_reference(frequency, temperatures, Path(outdir)/'Tcal',
                   title='Noise-diode reference | M{:02d} | {}'.format(beam,Path(Tcalfits).name))


def _plot_table(path, column, stem, unit, outdir, cal_bin=None):
    with fits.open(path) as hdus:
        values = hdus[2].data[column].copy()
        time = hdus[2].data['MJD'].copy()
        primary, sub = hdus[0].header, hdus[1].header
        freq = primary['OBSFREQ']-primary['OBSBW']/2 + np.arange(primary['OBSNCHAN'])*sub['CHAN_BW']
    if cal_bin is not None:
        values = values[:,cal_bin,:,:]
    common = shared_limits([values[:,i,:] for i in range(2)])
    for pol in range(2):
        waterfall(values[:,pol,:], time, freq, Path(outdir)/f'{stem}_p{pol}',
                  stem+' | '+Path(path).name, pol=pol, unit=unit, limits=common)


def plot_Tpsr(Tpsrfits, outdir='.'):
    _plot_table(Tpsrfits, 'TSYS_PSR', 'Tpsr2D', 'K', outdir)


def plot_off(calfits, outdir='.'):
    _plot_table(calfits, 'CAL', 'OFF2D', 'Readout (arb. units)', outdir, cal_bin=0)
