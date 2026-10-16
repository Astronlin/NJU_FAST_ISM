### Gridding with doppler tracking

import cygrid
import numpy as np
from read_crafts_spec import read_crafts_spec
from read_crafts_cal import read_crafts_cal
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord,SpectralCoord
from astropy.utils import iers
from FISH_utils import interp

#iers.conf.auto_download = False  
iers.conf.auto_max_age = None  

def gridding_with_doppler_tracking(speclist,center,imsize=[100,100],cell=1/60,rest_freq=1.420405751768,vmin=-700,vmax=700,dv=0.1,frame='lsrk',p=[0],doppler_tracking=True,outname='image.fits',psrflag=False,overwrite=False):

    c = 299792458.0/1e3 # km/s

    target_header = {
     'NAXIS':  2,
     'NAXIS1': imsize[0],
     'NAXIS2': imsize[1],
     'CTYPE1': 'RA---SIN',
     'CTYPE2': 'DEC--SIN',
     'CUNIT1': 'deg',
     'CUNIT2': 'deg',
     'CDELT1': -cell,
     'CDELT2': cell,
     'CRPIX1': imsize[0] / 2.,
     'CRPIX2': imsize[1] / 2.,
     'CRVAL1': center[0],
     'CRVAL2': center[1],
     }

    print(target_header)
    
    gridder       = cygrid.WcsGrid(target_header)
    FAST_beam     = 3/60  # unit: deg
    kernel_sigma  = FAST_beam/2/np.sqrt(8 * np.log(2))
    kernel_args   = ('gauss1d',(kernel_sigma,),3.*kernel_sigma,kernel_sigma/2.)
    gridder.set_kernel(*kernel_args)
    
    fmin   = (1-vmax/c)*rest_freq
    fmax   = (1-vmin/c)*rest_freq
    df     = dv/c*rest_freq          # unit: GHz;
    faxis = np.arange(fmin,fmax,df) # This is the frequency axis in the frame of the final cube
    
    speclist = np.loadtxt(speclist, dtype=str, ndmin=1)
    for specfits in speclist:
        try:
            if not psrflag:
                _,spec_hdu1 = read_crafts_spec(specfits)
            else:
                spec_hdu0,spec_hdu1,spec_hdu2 = read_crafts_cal(specfits)
        except:
            print('WARNING: BAD FITS FILE!! Please check '+ specfits)
        else:
            print('Gridding '+specfits)
            if not psrflag:
                spec_data   = np.nanmean(spec_hdu1.data['data'][:,:,p],axis=2) # p is a LIST 
                spec_fmin   = spec_hdu1.data['freq'][0]
                spec_chw    = spec_hdu1.data['chan_bw'][0]
                spec_nch    = spec_hdu1.data['nchan'][0]
                spec_fmax   = spec_fmin+spec_nch*spec_chw
                spec_freq   = np.arange(spec_fmin,spec_fmax,spec_chw)/1e3 # freq in topo
                spec_ra     = spec_hdu1.data['OBJ_RA']
                spec_dec    = spec_hdu1.data['OBJ_DEC']
                spec_mjd    = spec_hdu1.data['UTOBS']
            else:
                spec_data = np.nanmean(spec_hdu2.data['Tsys_psr'][:,p,:],axis=1) # p is a LIST
                obsfreq   = spec_hdu0.header['OBSFREQ']  # observing frequency
                obsbw     = spec_hdu0.header['OBSBW']    # observing band width
                spec_fmin = float(obsfreq - obsbw/2.)   # minimum frequency
                spec_fmax = float(obsfreq + obsbw/2.)   # maximum frequency
                spec_nch  = spec_hdu0.header['OBSNCHAN'] # observing channel number              
                spec_chw  = spec_hdu1.header['CHAN_BW']  # channel width
                spec_freq = np.arange(spec_fmin,spec_fmax,spec_chw)/1e3 # frequency array (GHz)
                spec_ra   = spec_hdu2.data['OBJ_RA']
                spec_dec  = spec_hdu2.data['OBJ_DEC']
                spec_mjd  = spec_hdu2.data['MJD']
                
            #####################################################################
            ####################### Doppler Tracking ############################
            # https://docs.astropy.org/en/stable/coordinates/spectralcoord.html #
            #####################################################################
            
            if doppler_tracking==True:
                
                print('Doing Doppler Tracking!')
                
                FAST = EarthLocation(x=-1668557.9202*u.m,y=5506865.964*u.m,z=2744946.0867*u.m)
                spec_time = Time(spec_mjd,format='mjd')
                
                if False:
                    data_to_regrid = np.tile(np.zeros_like(faxis,dtype=np.float32),(len(spec_ra),1))

                    # each spectrum (0.2 s)
                    for i in range(len(spec_ra)):

                        observer  = FAST.get_itrs(obstime=spec_time[i])
                        target    = SkyCoord(ra=spec_ra[i]*u.deg,dec=spec_dec[i]*u.deg,frame='icrs')
                        # There would be warning about the distance but it's okay for non-astrometry data.

                        sc_obs   = SpectralCoord(spec_freq*u.GHz,observer=observer,target=target)  
                        sc_frame = sc_obs.with_observer_stationary_relative_to(frame)

                        spec_freq_frame = sc_frame.value

                        ### Resample to the desired frequency axis ###
                        ### Updating the data_to_regrid array
                        data_to_regrid[i] = interp(spec_freq_frame,spec_data[i],faxis,kind='linear')
                else:
                    # each data file (0.2s x 2048)
                    observer  = FAST.get_itrs(obstime=Time(np.nanmean(spec_mjd),format='mjd'))
                    target    = SkyCoord(ra=np.nanmean(spec_ra)*u.deg,dec=np.nanmean(spec_dec)*u.deg,frame='icrs')
                    
                    sc_obs   = SpectralCoord(spec_freq*u.GHz,observer=observer,target=target)  
                    sc_frame = sc_obs.with_observer_stationary_relative_to(frame)
                    
                    spec_freq_frame = sc_frame.value
                    
                    data_to_regrid = interp(spec_freq_frame,spec_data,faxis,axis=1,kind='linear')
                    
                gridder.grid(spec_ra, spec_dec, data_to_regrid)
                
            else:
                print('No Doppler Tracking! Using the raw channel width.')
                freq_range  = (spec_freq>fmin) & (spec_freq<fmax) 
                spec_data   = np.nanmean(spec_data[:,freq_range,p],axis=2)
                spec_freq   = spec_freq[freq_range]
                gridder.grid(spec_ra, spec_dec, spec_data)

            ################################
            
    
    gridded_data = gridder.get_datacube()
    
    hdu = fits.PrimaryHDU(data=gridded_data)
    
    for key in target_header.keys():
        hdu.header[key] = target_header[key]
    
    hdu.header['CTYPE3']  = 'FREQ'
    hdu.header['CUNIT3']  = 'Hz'
    hdu.header['CDELT3']  = df*1e9
    hdu.header['CRPIX3']  = 0
    hdu.header['CRVAL3']  = fmin*1e9
    hdu.header['SPECSYS'] = frame.upper()
    hdu.header['RESTFRQ'] = rest_freq*1e9 
    hdu.header['BUNIT']   = 'K'

    hdu.writeto(outname,overwrite=overwrite)
