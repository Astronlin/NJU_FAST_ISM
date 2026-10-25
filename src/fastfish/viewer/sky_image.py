"""Bounded celestial FITS preview in the footprint's ICRS RA/Dec plane."""
import numpy as np

def load_sky_image(path,max_side=768):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.wcs.utils import pixel_to_skycoord,skycoord_to_pixel
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    with fits.open(path,memmap=True,do_not_scale_image_data=True) as hdus:
        selected=None
        for hdu in hdus:
            if hdu.header.get('NAXIS')==2 and hdu.header.get('NAXIS1',0)>1 and hdu.header.get('NAXIS2',0)>1:
                w=WCS(hdu.header)
                if w.has_celestial:selected=(hdu,w.celestial);break
        if selected is None:raise ValueError('Select a 2D FITS image with celestial WCS; cubes are not supported here')
        hdu,w=selected;ny,nx=hdu.shape
        t=np.linspace(0,1,129)
        xx=np.r_[t*(nx-1),np.full(129,nx-1),(1-t)*(nx-1),np.zeros(129)]
        yy=np.r_[np.zeros(129),t*(ny-1),np.full(129,ny-1),(1-t)*(ny-1)]
        sky=pixel_to_skycoord(xx,yy,w,origin=0).icrs
        ra=sky.ra.deg;dec=sky.dec.deg;finite=np.isfinite(ra)&np.isfinite(dec)
        if not finite.any():raise ValueError('Image WCS has no finite celestial boundary')
        x0,x1=float(ra[finite].min()),float(ra[finite].max());y0,y1=float(dec[finite].min()),float(dec[finite].max())
        if x1-x0>180:raise ValueError('Images crossing RA=0 require wrap-aware support; not supported yet')
        if x1<=x0 or y1<=y0:raise ValueError('Degenerate celestial image footprint')
        # Reproject at bounded display resolution; never create a full float copy of the source.
        width=min(max_side,nx);height=min(max_side,ny)
        dx=(x1-x0)/width;dy=(y1-y0)/height
        output=np.full((height,width),np.nan,dtype=np.float32)
        data=hdu.data;scale=hdu.header.get('BSCALE',1);offset=hdu.header.get('BZERO',0);blank=hdu.header.get('BLANK')
        for first in range(0,height,64):
            last=min(height,first+64)
            gx,gy=np.meshgrid(x0+(np.arange(width)+.5)*dx,y0+(np.arange(first,last)+.5)*dy)
            px,py=skycoord_to_pixel(SkyCoord(gx*u.deg,gy*u.deg,frame='icrs'),w,origin=0)
            good=np.isfinite(px)&np.isfinite(py)&(px>=0)&(py>=0)&(px<=nx-1)&(py<=ny-1)
            ix=np.floor(np.where(good,px,0)).astype(int);iy=np.floor(np.where(good,py,0)).astype(int)
            ix1=np.minimum(ix+1,nx-1);iy1=np.minimum(iy+1,ny-1)
            fx=np.where(good,px-ix,0);fy=np.where(good,py-iy,0)
            total=np.zeros(px.shape);weight=np.zeros(px.shape)
            for sx,sy,ww in ((ix,iy,(1-fx)*(1-fy)),(ix1,iy,fx*(1-fy)),(ix,iy1,(1-fx)*fy),(ix1,iy1,fx*fy)):
                values=np.asarray(data[sy,sx],dtype=float);valid=good&np.isfinite(values)
                if blank is not None:valid&=values!=blank
                values=values*scale+offset
                total+=np.where(valid,values,0)*ww;weight+=valid*ww
            chunk=output[first:last];np.divide(total,weight,out=chunk,where=weight>0)
        values=output[np.isfinite(output)]
        if not len(values):raise ValueError('Image reprojection contains no finite pixels')
        lo,hi=np.percentile(values,[1,99])
        if hi<=lo:hi=lo+1
        return dict(image=output,rect=(x0,y0,x1-x0,y1-y0),levels=(float(lo),float(hi)))
