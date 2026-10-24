## original version: from Weiwei Zhu email; Thu Nov 17;
## Modified version: Lingrui Lin; Nov 26, 2021
## read the CAL files of the FAST CRAFTS; can load in a list of fits files
## Modified version: Lingrui Lin; Nov 27, 2021
## save the ON-OFF spectrum object into a wcs1d-fits
## Modified version: Lingrui Lin; Jan 30, 2022
## simplify the code

# ------------ import necessary packages
import numpy as np
from astropy.io import fits
import sys
import gc
# ------------ Main Functions

def read_crafts_cal(calfile):
    '''
    Read CRAFTS High CAL files
    '''
    cal_hdu0,cal_hdu1,cal_hdu2 = fits.open(calfile)
    gc.collect()

    return cal_hdu0,cal_hdu1,cal_hdu2


# ---------- header information of the first fitsfile in fitslist
## comments: in the future, we need to check the consistency between headers of different observations
if __name__ == '__main__': 
    
    calfile = sys.argv[1]

    cal_hdu0,cal_hdu1,cal_hdu2 = read_crafts_cal(calfile)
    
    mjd  = cal_hdu2.data['mjd']
    data = cal_hdu2.data['cal']
    obsfreq  = cal_hdu0.header['OBSFREQ']                # observing frequency
    obsnchan = cal_hdu0.header['OBSNCHAN']               # observing channel number
    obsbw    = cal_hdu0.header['OBSBW']                  # observing band width
    fmin     = float(obsfreq - obsbw/2.)                # minimum frequency
    fmax     = float(obsfreq + obsbw/2.)                # maximum frequency
    nf       = obsnchan                                 # channel number
    df       = cal_hdu1.header['CHAN_BW']                # channel width
    
    # ----------- calculate the CAL_ON-OFF spectrum
    freq = np.arange(fmin,fmax,df)/1e3
    #datas     = data[:,:,:2,:].sum(axis=-1).sum(axis=-1)
    OFF       = np.average(data[:,0,:2,:].reshape(-1,nf),axis=0)
    ON        = np.average(data[:,1,:2,:].reshape(-1,nf),axis=0)
    NoiInj    = ON-OFF
    
    # ---------- save the CAL_ON-OFF spectrum
    #spec_NoiInj = Spectrum1D(flux=NoiInj*u.dimensionless_unscaled,spectral_axis=freq*u.MHz)
    #spec_NoiInj.write(filelist[:-5]+'_'+srcname+'_NoiInj.fits',overwrite=True)
    
    # ---------- smooth the CAL_ON-OFF spectrum for plot
    #spec_NoiInj_gsm  = gaussian_smooth(spec_NoiInj, stddev=10)
