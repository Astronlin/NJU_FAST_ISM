# transform the (x,y,z) coordinates of the feed to celestial frame

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from fastfish.io.read_crafts_spec import read_crafts_spec
from fastfish.io.read_crafts_cal import read_crafts_cal
import matplotlib.pyplot as plt
from astropy.utils import iers
from fastfish.calibration.interpolation import interp

#iers.conf.auto_download = False  
iers.conf.auto_max_age = None  

def reconstruct_pointing(spec_mjd, xyztable, beam=1):
    """Reconstruct ICRS pointing using time-dependent rotation in radians.

    FAST's feed-support reference specifies lengths in metres and angles in
    radians. Samples outside telemetry coverage receive NaN coordinates.
    """
    spec_mjd = np.asarray(spec_mjd, dtype=float)
    cabin_mjd = (Time(np.asarray(xyztable['SysTime'], dtype=str), format='iso') - 8*u.hr).mjd
    fields = [np.asarray(xyztable[name], dtype=float) for name in
              ['SDP_PhaPos_X', 'SDP_PhaPos_Y', 'SDP_PhaPos_Z', 'SDP_AngleM']]
    good = np.isfinite(cabin_mjd)
    for values in fields:
        good &= np.isfinite(values)
    cabin_mjd = cabin_mjd[good]
    fields = [values[good] for values in fields]
    order = np.argsort(cabin_mjd)
    cabin_mjd = cabin_mjd[order]
    fields = [values[order] for values in fields]
    if len(cabin_mjd) < 2 or np.any(np.diff(cabin_mjd) <= 0):
        raise ValueError('Pointing telemetry requires at least two unique, increasing timestamps')
    covered = np.isfinite(spec_mjd) & (spec_mjd >= cabin_mjd[0]) & (spec_mjd <= cabin_mjd[-1])
    ra = np.full(spec_mjd.shape, np.nan)
    dec = np.full(spec_mjd.shape, np.nan)
    if not covered.any():
        return ra, dec
    t = spec_mjd[covered]
    east, north, sky = [np.interp(t, cabin_mjd, values) for values in fields[:3]]
    angle = np.interp(t, cabin_mjd, np.unwrap(fields[3]))
    site = EarthLocation(x=-1668557.9202*u.m, y=5506865.964*u.m, z=2744946.0867*u.m)
    coord = SkyCoord(x=-north, y=-east, z=-sky, representation_type='cartesian',
                     frame='altaz', location=site, obstime=Time(t, format='mjd')).icrs
    beams = {'01':(0,0), '02':(5.74,0.00811), '03':(2.88,-4.97), '04':(-2.86,-4.98),
             '05':(-5.74,-0.0127), '06':(-2.88,4.97), '07':(2.86,4.98), '08':(11.5,0.0116),
             '09':(8.61,-4.96), '10':(5.75,-9.93), '11':(0.018,-9.94), '12':(-5.71,-9.96),
             '13':(-8.6,-4.99), '14':(-11.5,-0.03), '15':(-8.63,4.95), '16':(-5.77,9.93),
             '17':(-0.0181,9.94), '18':(5.73,9.95), '19':(8.61,4.98)}
    key = f'{int(beam):02d}'
    if key not in beams:
        raise ValueError('Beam must be between 1 and 19')
    dx, dy = beams[key]
    offset_east = dx*np.cos(angle) - dy*np.sin(angle)
    offset_north = dx*np.sin(angle) + dy*np.cos(angle)
    ra[covered] = (coord.ra.deg + offset_east/np.cos(coord.dec.radian)/60) % 360
    dec[covered] = coord.dec.deg + offset_north/60
    return ra, dec

def xyz2celestial(specfits,xyztable,beam=1,plot=False,psrflag=False):

    if not psrflag:

        spec_hdu0,spec_hdu1 = read_crafts_spec(specfits)

        spec_mjd  = spec_hdu1.data['utobs']
        spec_fmin = spec_hdu1.data['freq'][0]
        spec_chw  = spec_hdu1.data['chan_bw'][0]
        spec_nch  = spec_hdu1.data['nchan'][0]
        spec_fmax = spec_fmin+spec_nch*spec_chw
        spec_freq = np.arange(spec_fmin,spec_fmax,spec_chw)/1e3 # freq in topo
    else:
        spec_hdu0,spec_hdu1,spec_hdu2 = read_crafts_cal(specfits)
        spec_mjd = spec_hdu2.data['MJD']
    
    obj_ra, obj_dec = reconstruct_pointing(spec_mjd, xyztable, beam=beam)

    n = len(obj_ra)
    
    if not psrflag:
        spec_hdu1.data['OBJ_RA']  = obj_ra
        spec_hdu1.data['OBJ_DEC'] = obj_dec
        hduout = fits.HDUList([spec_hdu0,spec_hdu1])
    else:
        col_ra  = fits.Column(name='OBJ_RA',format="1D", array=obj_ra)
        col_dec = fits.Column(name='OBJ_DEC',format="1D", array=obj_dec)
        spec_hdu2.columns.add_col(col_ra)
        spec_hdu2.columns.add_col(col_dec)
        hduout = fits.HDUList([spec_hdu0,spec_hdu1,spec_hdu2])
    
    rootname  = specfits[:-5]
    hduout.writeto(rootname+'_icrs.fits', overwrite=True)
    hduout.close()

    if plot:
        fig = plt.figure()
        plt.plot(spec_hdu1.data['OBJ_RA'],spec_hdu1.data['OBJ_DEC'],color='black')
        plt.xlabel('R.A. [ICRS]')
        plt.ylabel('Dec. [ICRS]')
        plt.savefig(rootname+'_scan.pdf',format='pdf',bbox_inches='tight')
        plt.close()

def dp_tracking():

    ## This part is under construction.

    print('Star tracking.')