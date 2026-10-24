"""Viewport-dependent display sampling; never used to select science data."""
import numpy as np


def visible_indices(ra,dec,xrange,yrange,width,height,max_points=40000):
    """Return original indices, keeping one marker per occupied screen cell.

    Small visible selections retain every point. The overview is a coverage
    display, not a number-density estimator. Each occupied cell is represented.
    """
    ra=np.asarray(ra);dec=np.asarray(dec)
    lo,hi=sorted(xrange);bottom,top=sorted(yrange)
    ids=np.flatnonzero(np.isfinite(ra)&np.isfinite(dec)&(ra>=lo)&(ra<=hi)&(dec>=bottom)&(dec<=top))
    if len(ids)<=max_points:return ids
    aspect=max(1.,float(width))/max(1.,float(height));nx=max(1,min(max_points,int(np.sqrt(max_points*aspect))));ny=max(1,max_points//nx)
    x=np.minimum(nx-1,((ra[ids]-lo)/max(hi-lo,np.finfo(float).tiny)*nx).astype(np.int64))
    y=np.minimum(ny-1,((dec[ids]-bottom)/max(top-bottom,np.finfo(float).tiny)*ny).astype(np.int64))
    _,first=np.unique(y*nx+x,return_index=True)
    return ids[np.sort(first)]
