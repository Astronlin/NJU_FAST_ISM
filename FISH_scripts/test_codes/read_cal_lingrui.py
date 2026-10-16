## original version: from Weiwei Zhu email; Thu Nov 17;
## Modified version: Lingrui Lin; Nov 26, 2021
## read the CAL files of the FAST CRAFTS; can load in a list of fits files
## Modified version: Lingrui Lin; Nov 27, 2021
## save the ON-OFF spectrum object into a wcs1d-fits

# ------------ import necessary packages
import numpy as np
import sys
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from specutils import Spectrum1D
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)
from pylab import *
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
# ------------ set the multiprocessing parameters
#ncores = mp.cpu_count()
#nuse   = 5
#print('Total cpu cores:',ncores)
#print('Used cores:',nuse)
#def readhdu(fitsfile):
#    hdu = fits.open(fitslist[i])
#    hdu0.append(hdu[0])
#    hdu1.append(hdu[1])
#    hdu2.append(hdu[2])
#    mjds.append(hdu[2].data['MJD'])
#    data.append(hdu[2].data['CAL']) 
# ------------ open the CAL file and read in the HDUs
filelist = sys.argv[1]
srcname  = sys.argv[2]
fitslist = np.loadtxt(filelist,dtype='str')
fitslist = [selected for selected in fitslist if srcname in selected]
#os.system('mv *M{}* M{}'.format(str(i).rjust(2,'0'),str(i).rjust(2,'0')))
hdu0 = []
hdu1 = []
hdu2 = []
mjds = []
data = []
for i in range(len(fitslist)):
    #print('Reading {}'.format(fitslist[i]))
    hdu = fits.open(fitslist[i])
    hdu0.append(hdu[0])
    hdu1.append(hdu[1])
    hdu2.append(hdu[2])
    mjds.append(hdu[2].data['MJD'])                    # MJD
    data.append(hdu[2].data['CAL'])                    # CAL info

data = np.concatenate(data,axis=0)
mjds = np.concatenate(mjds,axis=0)

# ---------- header information of the first fitsfile in fitslist
## comments: in the future, we need to check the consistency between headers of different observations
obsfreq  = hdu0[0].header['OBSFREQ']                # observing frequency
obsnchan = hdu0[0].header['OBSNCHAN']               # observing channel number
obsbw    = hdu0[0].header['OBSBW']                  # observing band width
fmin     = float(obsfreq - obsbw/2.)                # minimum frequency
fmax     = float(obsfreq + obsbw/2.)                # maximum frequency
nf       = obsnchan                                 # channel number
df       = hdu1[0].header['CHAN_BW']                # channel width

# ----------- calculate the CAL_ON-OFF spectrum
data_freq = np.arange(fmin,fmax,df)/1e3
#datas     = data[:,:,:2,:].sum(axis=-1).sum(axis=-1)
OFF       = np.average(data[:,0,:2,:].reshape(-1,nf),axis=0)
ON        = np.average(data[:,1,:2,:].reshape(-1,nf),axis=0)
NoiInj    = ON-OFF

# ---------- save the CAL_ON-OFF spectrum
spec_NoiInj = Spectrum1D(flux=NoiInj*u.dimensionless_unscaled,spectral_axis=data_freq*u.MHz)
spec_NoiInj.write(filelist[:-5]+'_'+srcname+'_NoiInj.fits',overwrite=True)

# ---------- smooth the CAL_ON-OFF spectrum for plot
spec_NoiInj_gsm  = gaussian_smooth(spec_NoiInj, stddev=10)

# ------------- plot the results
fig = plt.figure()
gs = fig.add_gridspec(2,hspace=0,height_ratios=np.array([1.5,1]))
ax = gs.subplots(sharex=True,sharey=False)
ax[0].step(data_freq,ON,where='mid',color='red',label='CAL_ON')
ax[0].step(data_freq,OFF,where='mid',color='blue',label='CAL_OFF')
ax[1].step(data_freq,NoiInj,where='mid',color='green',label='CAL_ON-CAL_OFF')
ax[1].step(spec_NoiInj_gsm.spectral_axis,spec_NoiInj_gsm.flux,where='mid',color='orange',label='CAL_ON-CAL_OFF smooth')
ax[0].legend()
ax[1].legend()
ax[0].set_ylim(0,20)
ax[1].set_xlabel('Frequency [GHz]',fontweight='bold')
ax[0].set_title(filelist[:-5]+'_'+srcname)
outfig = filelist[:-5]+'_'+srcname+'_CAL.pdf'
plt.savefig(outfig,bbox_inches='tight',format='pdf')
plt.close()
#os.system('xdg-open {}'.format(outfig))

'''
data = npzfile['data']
print(data.shape)

datainfo = npzfile['dimensions']
print(datainfo)

header = npzfile['header']
print(header)

'''


#tsamp    = hdu1[0].header['TBIN']                      # the sampling time
#nsubint  = hdu1[0].header['NAXIS2']                  # number of subintegrations
#samppersubint = int(hdu1[0].header['NSBLK'])        # sampling per subintegration
#nsamp = nsubint * samppersubint                  # total sampling number
#sourename = hdu0[0].header['SRC_NAME']              # sourcename
#nbits = hdu0[0].header['BITPIX']                    # bits of datum
#header = dict(hdu0[0].header + hdu1[0].header)         # combination of the first two headers
#duration = tsamp * nsamp                         # time duration of the observation

# -------------- print necessary information
#print('obsfreq: {}MHz'.format(obsfreq))
#print('obsnchan: {}'.format(obsnchan))
#print('tsamp: {}s'.format(tsamp))
#print('duration: {}s'.format(duration,'s'))
#print(mjds)
#print(data.shape)
# data.shape = (64, 2, 4, 4096) :            64: subintegrations, 2: ON+OFF, 4: polarization, 4096: channel number
#print hdu1.data['TSUBINT']
#print hdu1.data['OFFS_SUB']
#sys.exit(0)
