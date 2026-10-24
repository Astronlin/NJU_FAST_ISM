"""Vectorized per-spectrum frequency-frame factors, without spectral resampling."""
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation,SkyCoord,SpectralCoord,CartesianDifferential,CartesianRepresentation,frame_transform_graph
from astropy.time import Time
SITE=EarthLocation(x=-1668557.9202*u.m,y=5506865.964*u.m,z=2744946.0867*u.m)


def frequency_scales(mjd,ra_deg,dec_deg,frame='lsrk',rest_freq_ghz=1.420405751768,mode='spectrum'):
    t,ra,dec=[np.asarray(x,dtype=float) for x in (mjd,ra_deg,dec_deg)]
    if t.ndim != 1 or t.shape != ra.shape or t.shape != dec.shape or not len(t) or not all(np.isfinite(x).all() for x in (t,ra,dec)):
        raise ValueError('Finite matching one-dimensional times and directions required')
    if mode == 'file':
        target=SkyCoord(float(ra.mean())*u.deg,float(dec.mean())*u.deg,frame='icrs')
        observer=SITE.get_itrs(obstime=Time(float(t.mean()),format='mjd'))
        value=SpectralCoord(rest_freq_ghz*u.GHz,observer=observer,target=target).with_observer_stationary_relative_to(frame).value
        return np.full(len(t),float(value)/rest_freq_ghz)
    if mode != 'spectrum':raise ValueError('Doppler mode must be spectrum or file')
    # Explicit array-shaped zero velocities avoid Astropy 5's scalar
    # attach_zero_velocities fallback, and match its scalar infinite-distance
    # convention. Zero ITRS velocity means a site fixed to the rotating Earth.
    zeros=CartesianDifferential(np.zeros((3,len(t)))*u.km/u.s)
    target=SkyCoord(ra*u.deg,dec*u.deg,distance=np.full(len(t),1e6)*u.kpc,frame='icrs').frame
    target=target.realize_frame(target.cartesian.with_differentials(zeros))
    observer=SITE.get_itrs(obstime=Time(t,format='mjd'))
    observer=observer.realize_frame(observer.cartesian.with_differentials(zeros))
    frequency=np.full(len(t),rest_freq_ghz)*u.GHz
    reference_cls=frame_transform_graph.lookup_name(frame)
    if reference_cls is None:raise ValueError('Unknown stationary frame: '+str(frame))
    stationary=reference_cls(CartesianRepresentation(np.zeros((3,len(t)))*u.km,differentials=zeros))
    value=SpectralCoord(frequency,observer=observer,target=target).with_observer_stationary_relative_to(stationary).value
    return np.asarray(value)/rest_freq_ghz
