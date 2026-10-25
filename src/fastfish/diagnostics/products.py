"""Compact cube diagnostics; display averages do not change science products."""
from pathlib import Path
import json
import warnings
import numpy as np
from astropy.io import fits


def plot_products(output_dir, name):
    import matplotlib.pyplot as plt
    root=Path(output_dir)
    with fits.open(root/(name+'_cube.fits'),memmap=True) as h:
        data=h[0].data; hdr=h[0].header
        f=hdr['CRVAL3']+np.arange(len(data))*hdr['CDELT3']
        velocity=(1-f/hdr['RESTFRQ'])*299792.458
        centre=data[np.argmin(abs(velocity))].copy()
        # The wing window is diagnostic, not assumed to be emission-free.
        ids=np.flatnonzero((abs(velocity)>=80)&(abs(velocity)<=100))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            wing=np.nanmean(data[ids],axis=0)
            profile=np.array([np.nanmean(data[i]) for i in range(len(data))])
    with fits.open(root/'coverage.fits',memmap=True) as h:
        coverage=h[0].data.copy()
    with fits.open(root/'rfi_contribution_fraction.fits',memmap=True) as h:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            marks=np.array([np.nanmean(h[0].data[i]) for i in range(len(velocity))])
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ax,image,title,unit in zip(axes.flat[:3],[centre,wing,coverage],
        ['Near 0 km/s','Mean |v|=80–100 km/s; may contain sky/RFI','Geometric coverage'],['K','K','weight']):
        im=ax.imshow(image,origin='lower',cmap='viridis'); ax.set_title(title)
        ax.set_xlabel('RA pixel'); ax.set_ylabel('Dec pixel'); fig.colorbar(im,ax=ax,label=unit)
    axes[1,1].plot(velocity,profile,lw=1,label='Spatial mean')
    axes[1,1].set(xlabel='LSRK radio velocity (km/s)',ylabel='K')
    ax2=axes[1,1].twinx(); ax2.plot(velocity,marks,color='red',alpha=.6,label='Mean marked fraction')
    ax2.set_ylabel('RFI/candidate response fraction',color='red'); ax2.set_ylim(0,1)
    fig.suptitle(name+' | retained RFI values; baseline inferred inside protected H I')
    dest=root/'diagnostics.png'; fig.savefig(dest,dpi=160); plt.close(fig)
    return dest
