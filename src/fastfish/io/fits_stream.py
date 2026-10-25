"""Read fixed-width FAST binary tables without eager FITS column conversion.

Astropy 5 may construct logical columns while materializing a whole table.
For FAST's megabyte-wide rows that strided access provokes enormous read-ahead,
even when asking only for times. Header parsing remains with Astropy; numeric
table fields are read as their explicit big-endian FITS dtype. No scaled or
variable-length numeric columns are silently accepted.
"""
from contextlib import contextmanager
from pathlib import Path
import mmap
import os
import numpy as np
from astropy.io import fits


@contextmanager
def fixed_table(path, extension=1, random_access=False, sequential_access=False):
    """Yield (numeric table view, HDU headers); views expire on context exit.

    Logical columns remain raw FITS bytes and must not be interpreted through
    this helper. The pipeline uses only unscaled fixed-width numeric columns.
    """
    if random_access and sequential_access:raise ValueError('Choose one FITS access hint')
    with fits.open(path,memmap=True) as h:
        hd=h[extension]
        for i,column in enumerate(hd.columns,1):
            if 'P' in str(column.format) or 'Q' in str(column.format):
                raise ValueError('Variable-length FITS fields are not supported')
            if hd.header.get('TSCAL%d'%i,1)!=1 or hd.header.get('TZERO%d'%i,0)!=0:
                raise ValueError('Scaled FITS fields need explicit conversion')
        dtype=hd.columns.dtype.newbyteorder('>')
        if dtype.itemsize!=int(hd.header['NAXIS1']):raise ValueError('FITS row length differs from dtype')
        with open(path,'rb',buffering=0) as f:
            with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
                if hasattr(mm,'madvise'):
                    if random_access:mm.madvise(mmap.MADV_RANDOM)
                    elif sequential_access:mm.madvise(mmap.MADV_SEQUENTIAL)
                table=np.ndarray((int(hd.header['NAXIS2']),),dtype=dtype,buffer=mm,offset=hd._data_offset)
                yield table,h
                del table


def spectral_metadata(path, read_times=False, time_reference='start'):
    """Read scalar prefixes; all effective science times refer to mid-integration.

    UTOBS is preserved as raw metadata. For start-tagged data, shift by each
    row's EXPOSURE/2. Already midpoint-tagged data can opt out explicitly.
    """
    if time_reference not in ('start','midpoint'):raise ValueError('Science time_reference must be start or midpoint')
    with fits.open(path,memmap=True) as h:
        hd=h[1];dtype=hd.columns.dtype;rowbytes=int(hd.header['NAXIS1'])
        n=int(hd.header['NAXIS2']);offset=hd._data_offset;fields=dtype.fields
        if n<1:raise ValueError('Empty science table')
        prefix_length=int(fields['DATA'][1]);time_dtype,time_offset=fields['UTOBS'][:2]
        if time_reference=='start' and 'EXPOSURE' not in fields:
            raise ValueError('Start-tagged science data require EXPOSURE to calculate integration midpoints')
        with open(path,'rb',buffering=0) as f:
            if hasattr(os,'posix_fadvise'):os.posix_fadvise(f.fileno(),0,0,os.POSIX_FADV_RANDOM)
            first=os.pread(f.fileno(),prefix_length,offset)
            last=os.pread(f.fileno(),prefix_length,offset+(n-1)*rowbytes)
            def scalar(name,buffer=first):
                dt,position=fields[name][:2]
                return np.ndarray((),dtype=dt.newbyteorder('>'),buffer=buffer,offset=position).item()
            raw=np.asarray([float(scalar('UTOBS')),float(scalar('UTOBS',last))])
            exposure=None
            if 'EXPOSURE' in fields:
                exposure=np.asarray([float(scalar('EXPOSURE')),float(scalar('EXPOSURE',last))])
                if not np.isfinite(exposure).all() or (exposure<=0).any():raise ValueError('Science EXPOSURE must be finite and positive')
            effective=raw+exposure/(2*86400) if time_reference=='start' else raw
            info={'path':str(path),'bytes':Path(path).stat().st_size,'rows':n,'mjd_first_last':effective.tolist(),
                  'raw_mjd_first_last':raw.tolist(),'raw_timestamp_reference':time_reference,'timestamp_reference':'midpoint',
                  'exposure_first_last_s':exposure.tolist() if exposure is not None else None,
                  'freq_first_ghz':float(scalar('FREQ'))/1e3,'native_df_hz':float(scalar('CHAN_BW'))*1e6,
                  'channels':int(scalar('NCHAN'))}
            if read_times:
                with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
                    if hasattr(mm,'madvise'):mm.madvise(mmap.MADV_RANDOM)
                    times=np.ndarray((n,),dtype=time_dtype.newbyteorder('>'),buffer=mm,offset=offset+time_offset,strides=(rowbytes,)).astype(float)
                    if time_reference=='start':
                        dt,position=fields['EXPOSURE'][:2]
                        durations=np.ndarray((n,),dtype=dt.newbyteorder('>'),buffer=mm,offset=offset+position,strides=(rowbytes,)).astype(float)
                        if not np.isfinite(durations).all() or (durations<=0).any():raise ValueError('Science EXPOSURE must be finite and positive')
                        times+=durations/(2*86400)
                return info,times
    return info


def contiguous_rows(ids):
    """Only consecutive row selections qualify for streaming read-ahead."""
    return len(ids)>1 and bool(np.all(np.diff(ids)==1))
