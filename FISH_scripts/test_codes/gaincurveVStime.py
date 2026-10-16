## To check whether the gaincurve of the Tcal injection varies with time
## Lingrui, Dec 16, 2021

# ------------------- import necessary packages
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

data    = np.concatenate(data,axis=0)
mjds    = np.concatenate(mjds,axis=0)
nsubint = len(mjds)
# ---------- header information of the first fitsfile in fitslist
## comments: in the future, we need to check the consistency between headers of different observations

obsfreq  = hdu0[0].header['OBSFREQ']                # observing frequency
obsnchan = hdu0[0].header['OBSNCHAN']               # observing channel number
obsbw    = hdu0[0].header['OBSBW']                  # observing band width
fmin     = float(obsfreq - obsbw/2.)                # minimum frequency
fmax     = float(obsfreq + obsbw/2.)                # maximum frequency
nf       = obsnchan                                 # channel number
df       = hdu1[0].header['CHAN_BW']                # channel width

data_freq = np.arange(fmin,fmax,df)/1e3

# ----------- calculate the CAL_ON-OFF spectrum (do not average the subintegrations)
ON     = data[:,1,:2,:]
OFF    = data[:,0,:2,:]
NoiInj = np.sum(ON-OFF,axis=1)

# ------------ plot the gain curves
fig = plt.figure()
gs  = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
ax00   = fig.add_subplot(gs[0,0])
ax01   = fig.add_subplot(gs[0,1],sharey=ax00)
ax10   = fig.add_subplot(gs[1,0],sharex=ax00)
Noimap = ax00.imshow(NoiInj,cmap='jet',aspect='auto',vmin=0,vmax=1)
ax10.step(np.arange(0,nf),np.average(NoiInj,axis=0),where='mid')
ax01.step(np.average(NoiInj,axis=1),np.arange(0,nsubint),where='mid')
ax10.set_ylim(0,1.3)
ax00.tick_params(labelbottom=False)
ax01.tick_params(labelleft=False)
ax00.set_ylabel('Subintegration')
ax01.set_xlabel('ON-OFF')
ax10.set_xlabel('Channel')#,fontweight='bold')
ax10.set_ylabel('ON-OFF')

#plt.colorbar(mappable=Noimap,pad=0.0,ax=ax00,location='top',aspect=25)
plt.savefig('gaincurveVStime.pdf')
plt.close()

# -------------- History codes
##ax  = gs.subplots(sharex=True,sharey=False)
##
##Noimap = ax[0][0].imshow(NoiInj,cmap='jet',aspect='auto',vmin=0,vmax=1)
##ax[1][0].step(np.arange(0,nf),np.average(NoiInj,axis=0),where='mid')
##ax[0][0].set_ylabel('Subintegration')
##ax[1][0].set_xlabel('Channel')#,fontweight='bold')
##ax[1][0].set_ylabel('ON-OFF')
##ax[1][0].set_ylim(0,1.5)


