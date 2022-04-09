## Extract the CAL_ON, CAL_OFF from the PSR data
## Lastest version: Mar 28, 2022

import numpy as np
from astropy.io import fits
from astropy.time import Time as ATime
from datetime import timedelta
import sys
#import time

def extract_crafts_cal(filename):
    
    print('Extracting CAL from '+filename)
    
    hdulist = fits.open(filename)
    hdu0    = hdulist[0]
    hdu1    = hdulist[1]

    #obsfreq       = hdu0.header['OBSFREQ']          # center frequency
    #obsnchan      = hdu0.header['OBSNCHAN']         # number of frequency channels
    #obsbw         = hdu0.header['OBSBW']            # total bandwidth
    #sourename     = hdu0.header['SRC_NAME']
    #nbits         = hdu0.header['BITPIX']
    #nf            = obsnchan
    #df            = hdu1.header['CHAN_BW']          # channel width
    tsamp         = hdu1.header['TBIN']             # sampling time
    nsubint       = hdu1.header['NAXIS2']           # number of subintegrations
    samppersubint = int(hdu1.header['NSBLK'])       # sampling numbers per subint
    nsamp         = nsubint * samppersubint         # total samplings
    #header        = dict(hdu0.header + hdu1.header)
    duration      = tsamp * nsamp                   # total observing time (s)
    Tsubint       = samppersubint * tsamp           # 0.100663296 (s)

    Nbin          = 2                          # number of bins to average along subint
    Tconst        = Nbin * Tsubint             # (28672+20480)/2*4*1024*2*1e-9 s
    samppersubint = int(Tconst / tsamp)        # new sampling number per subint
    Nsub          = int(nsamp / samppersubint) # new subint number
    Tsubint       = Tconst                     # new tsubint
    
    obsmjd = int(hdu0.header['STT_IMJD'])     # START INT MJD (UTC days)
    tstart = float(hdu0.header['STT_SMJD']) \
            + float(hdu0.header['STT_OFFS'])  # START time + START time offset (s)
    mjd    = obsmjd + tstart/86400            # FLOAT MJD at the start of the observation
    
    obsdate = ATime(obsmjd, format='mjd').datetime
    obstime = timedelta( seconds = float(tstart) )
    obsutc  = obsdate + obstime
    
    #header['Tsubint']  = Tsubint
    #header['Duration'] = duration
    #header['Obsmjd']   = str(mjd)
    #header['obsutc']   = obsutc
    
    data0     = np.copy(hdu1.data['data'][:,:,:,:,:])
    _,_,p,n,_ = data0.shape
    
    data  = data0.reshape(Nsub,-1,Nbin,p,n).mean(axis=1, dtype=np.float32) # Caveat: data type of AB and BA?
                                                                          # group the 128 subintegrations two-by-two;
                                                                          # group the 2*1024 sampling two-by-two;
                                                                          # average in each subintegration group
    Tarray = mjd + (np.arange(Nsub) + 0.5) * Tconst / 86400 # center time (float day) of each subint
    
    #delnames = ['DATA', 'DAT_OFFS', 'DAT_SCL', 'DAT_WTS', 'DAT_FREQ']
    #datadict = {}
    #for name in hdu1.data.names:
    #    if name not in delnames:
    #        datadict[name] = hdu1.data[name]
    
    hdu1.columns.del_col('DATA')
    hdu1.columns.del_col('DAT_OFFS')
    hdu1.columns.del_col('DAT_SCL')
    hdu1.columns.del_col('DAT_WTS')
    hdu1.columns.del_col('DAT_FREQ')
    fits.BinTableHDU.update(hdu1)
    #for name in datadict:
    #    hdu1.data[name] = datadict[name]
    
    CAL    = fits.Column(name='CAL', format="%dE" % (Nbin*p*n), array=data, dim=str((n, p, Nbin))) # OFF, ON, OFF, ON...
    TIME   = fits.Column(name='MJD', format="D", array=Tarray)
    cols   = fits.ColDefs([TIME, CAL])
    caltab = fits.BinTableHDU.from_columns(cols)
    hdulist.append(caltab)
    
    rootname = filename[:-5].split('/')[-1]
    decstr   = rootname.split('-M')[0]
    beamnum  = rootname.split('-M')[-1]
    newname  = 'HighCAL/'+ decstr + '-M' + beamnum + "_" + obsutc.strftime("%y%m%d") + "_" + obsutc.strftime("%H%M%S")
    hdulist.writeto(newname + '_cal.fits', overwrite=True)
    print("CAL SUCCESSFULLY EXTRACTED! "+ newname + '_cal.fits')

if __name__ == '__main__':
    filename = sys.argv[1]
    extract_crafts_cal(filename)
