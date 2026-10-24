# read in the Spec-backend data
# Lingrui

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import gc

def read_crafts_spec(specfile):
    '''
    Read Spec Files
    '''
    spec_hdu0,spec_hdu1 = fits.open(specfile)
    gc.collect()
    return spec_hdu0,spec_hdu1


if __name__ == '__main__':

    specfile = sys.argv[1]
    
    spec_hdu0,spec_hdu1 = read_crafts_spec(specfile)
    
    utobs = spec_hdu1.data['utobs']
    f0    = spec_hdu1.data['freq'][0]
    nchan = spec_hdu1.data['nchan'][0]
    chanw = spec_hdu1.data['CHAN_BW'][0]
    data  = spec_hdu1.data['data']
    fmin  = f0
    fmax  = f0+nchan*chanw
    freq  = np.arange(fmin,fmax,chanw)
    
    data_p0 = data[:,:,0] # four polarizations
    data_p1 = data[:,:,1]
    data_p2 = data[:,:,2]
    data_p3 = data[:,:,3]
    
    # ----- test plot
    plt.figure()
    plt.plot(freq,np.average(data_p0,axis=0),label='p0')
    plt.plot(freq,np.average(data_p1,axis=0),label='p1')
    plt.plot(freq,np.average(data_p2,axis=0),label='p2')
    plt.plot(freq,np.average(data_p3,axis=0),label='p3')
    plt.legend()
    plt.savefig(specfile+'.pdf',format='pdf')
    plt.close()
    
