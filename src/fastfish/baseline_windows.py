"""Validated LSRK baseline exclusions, shared by preview and reduction recipes."""
import numpy as np

def normalize_windows(windows):
    if len(windows)>64:raise ValueError('At most 64 extra baseline exclusion intervals')
    result=[]
    for pair in windows:
        if len(pair)!=2:raise ValueError('Each exclusion needs two velocity limits')
        lo,hi=map(float,pair)
        if not np.isfinite([lo,hi]).all() or lo>=hi:raise ValueError('Exclusions require finite increasing velocity limits')
        result.append([lo,hi])
    merged=[]
    for lo,hi in sorted(result):
        if merged and lo<=merged[-1][1]:merged[-1][1]=max(hi,merged[-1][1])
        else:merged.append([lo,hi])
    return merged

def parse_windows(text):
    rows=[row.strip() for row in text.replace(';','\n').splitlines() if row.strip()]
    try:return normalize_windows([row.replace(':',',').split(',') for row in rows])
    except (TypeError,ValueError) as e:raise ValueError('Extra exclusions: enter one lo,hi pair per line (km/s). '+str(e)) from e

def format_windows(windows):
    return '\n'.join('{:.8g},{:.8g}'.format(lo,hi) for lo,hi in normalize_windows(windows))
