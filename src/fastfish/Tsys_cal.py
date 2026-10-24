# Generate the Tsys_psr

import numpy as np
import sys
from .read_crafts_cal import read_crafts_cal
from .read_crafts_spec import read_crafts_spec
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time as ATime
from scipy.interpolate import interp1d
from scipy import ndimage
from .FISH_utils import interp

def gaincal(callist,outname='Tsys_psr',do_cal=True,Tcalfits='Tcal/20201014/CAL.20201014.low.W.fits',beam=1,sm_nf=16,sm_nt=300,plot=False):
    
    callist = np.loadtxt(callist,dtype=str)
    
    cal_mjd = []
    cal_on  = []
    cal_off = []
    
    for calfile in callist:
        try:
            hdu0,hdu1,hdu2 = read_crafts_cal(calfile)
            cal_mjd.append(hdu2.data['mjd'])
            cal_on.append(hdu2.data['cal'][:,1,:,:])
            cal_off.append(hdu2.data['cal'][:,0,:,:])
        except:
            print('WARNING: BAD FITS FILE!! Please check '+ calfile)
    
    cal_mjd = np.concatenate(cal_mjd,axis=0)
    cal_on  = np.concatenate(cal_on,axis=0)
    cal_off = np.concatenate(cal_off,axis=0)
    onoff   = (cal_on-cal_off)/20480*(20480+28672)/2
    nsubint = len(cal_mjd)
    
    del cal_on

    obsfreq  = hdu0.header['OBSFREQ']                # observing frequency
    obsnchan = hdu0.header['OBSNCHAN']               # observing channel number
    obsbw    = hdu0.header['OBSBW']                  # observing band width
    fmin     = float(obsfreq - obsbw/2.)             # minimum frequency
    fmax     = float(obsfreq + obsbw/2.)             # maximum frequency
    nf       = obsnchan                              # channel number
    df       = hdu1.header['CHAN_BW']                # channel width
    cal_freq = np.arange(fmin,fmax,df)/1e3           # frequency array; MHz to GHz
    p        = onoff.shape[1]                        # n_polarization
    
    onoff_interp      = interp(cal_mjd    , onoff      , cal_mjd  , kind = 'linear' , axis = 0 , smooth = sm_nt , smoothmode = 'nearest') # smooth along mjd
    
    if do_cal:
        
        _,noisehdu1 = read_crafts_spec(Tcalfits)
        noise_data  = np.transpose(noisehdu1.data['tcal'][0,beam-1,:,:],axes=(1,0))
        print('There are {} polarizations in the Tcal noise.'.format(noise_data.shape[0]))
        if noise_data.shape[0]<p: 
        # add NaN array when the length of the polarization axis of the noise data is smaller than 4.
            noise_data = np.concatenate((noise_data,np.nan*np.ones(shape=(p-noise_data.shape[0],noise_data.shape[1]))),axis=0)
        noise_freq  = noisehdu1.data['freq'][0]/1e3 # MHz to GHz

        noise_data_interp = interp(noise_freq , noise_data , cal_freq , kind = 'linear' , axis = 1 , smooth = sm_nf  , smoothmode = 'nearest') # smooth along frequency

        Tpsr              = noise_data_interp/onoff_interp*cal_off ### Tsys_psr

        TSYS_PSR = fits.Column(name='TSYS_PSR', format="%dE" % (p*nf), array=Tpsr,dim=str((nf, p)))
        TIME     = fits.Column(name='MJD', format="D", array=cal_mjd)
        cols     = fits.ColDefs([TSYS_PSR,TIME])
        tsystab  = fits.BinTableHDU.from_columns(cols)
        hduout   = fits.HDUList([hdu0,hdu1,tsystab])

        cal_date = ATime(cal_mjd,format='mjd').datetime
        newname = outname + '-M{}'.format(str(beam).rjust(2,'0')) + '_' + cal_date[0].strftime("%y%m%d") + "_" + cal_date[0].strftime("%H%M%S") + 'to' \
                  + cal_date[-1].strftime("%y%m%d") + "_" + cal_date[-1].strftime("%H%M%S")
        hduout.writeto(newname+'.fits', overwrite=True)
        hduout.close()

    if plot:
        from .FISH_diagnostics import waterfall, shared_limits, tcal_reference
        from pathlib import Path
        context = '{} | M{:02d}'.format(Path(calfile).stem.split('-M')[0], beam)
        common = shared_limits([onoff[:, i, :] for i in range(2)])
        for i in range(2):
            waterfall(onoff[:, i, :], cal_mjd, cal_freq*1e3,
                      'onoff2D_p{}'.format(i), 'Diode ON−OFF | '+context,
                      pol=i, unit='Readout (arb. units)', reference=onoff_interp[:, i, :],
                      reference_label='Median after {}-sample channelwise smoothing'.format(sm_nt), limits=common,
                      note='Duty-cycle corrected; no RFI mask applied.')
        if do_cal:
            tcal_reference(noise_freq*1e3, noise_data[:2], 'Tcal_resample',
                           target_freq_mhz=cal_freq*1e3, processed=noise_data_interp[:2],
                           title='Tcal resampling | '+context)
            common = shared_limits([Tpsr[:, i, :] for i in range(2)])
            for i in range(2):
                waterfall(Tpsr[:, i, :], cal_mjd, cal_freq*1e3,
                          'Tpsr2D_p{}'.format(i), 'PSR temperature reference | '+context,
                          pol=i, unit='K', limits=common,
                          note='Full-band overview; invalid values are gray.')
            hi_limits = np.array([-700,700])/3e5*1.420405751768 + 1.420405751768
            hi = (cal_freq>hi_limits[0]) & (cal_freq<hi_limits[1])
            if np.count_nonzero(hi)>1:
                zoom_limits = shared_limits([Tpsr[:, i, :][:, hi] for i in range(2)])
                for i in range(2):
                    waterfall(Tpsr[:, i, :][:, hi], cal_mjd, cal_freq[hi]*1e3,
                              'Tpsr_HI_p{}'.format(i), 'PSR temperature reference — H I zoom | '+context,
                              pol=i, unit='K', limits=zoom_limits,
                              note='PSR channel spacing is not the spectral-backend resolution.')


    return  

if __name__ == '__main__':
    callist = sys.argv[1]
    plot = sys.argv[2]
    gaincal(callist,plot)
