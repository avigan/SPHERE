# -*- coding: utf-8 -*-

'''
Aperture/pupil utility functions for pyZELDA
'''

import numpy as np
import collections


def coordinates(dim, size, diameter=False, strict=False, center=(), cpix=False, normalized=True, outside=np.nan, polar=True):
    '''
    Returns rho,theta coordinates defining a circular aperture
    
    Parameters
    ----------
    dim : int
        Size of the output array
    
    size : int
        Size of the disk, representing either the radius (default) or the diameter
    
    diameter : bool, optional
        Specify if input size is the diameter. Default is 'False'
    
    strict : bool optional
        If set to Trye, size must be strictly less than (<), instead of less
        or equal (<=). Default is 'False'
    
    center : sequence, optional
        Specify the center of the disc. Default is '()', i.e. the center of the array
    
    cpix : bool optional
        If set to True, the disc is centered on pixel at position (dim//2,dim//2).
        Default is 'False', i.e. the disc is centered between 4 pixels
    
    normalized : bool optional
        Determines if the rho coordinates are normalized to 1. Default is True
    
    outside : float, optional
        Value used to fill the array outside of the aperture. Default is np.nan
    
    polar : bool, optional
        Return polar coordinates. Default is True
    
    Returns
    -------
    rho, theta : arrays
        Arrays containing the rho and theta coordinates of the aperture
    '''
    
    if (diameter is True):
        rad = size/2
    else:
        rad = size

    if (len(center) == 0):
        if (cpix is False):
            cx = (dim-1) / 2
            cy = (dim-1) / 2
        else:
            cx = dim // 2
            cy = dim // 2
    elif (len(center) == 2):
        cx = center[0]
        cy = center[1]
    else:
        raise ValueError('Error, you must pass 2 values for center')
        return None

    x = np.arange(dim, dtype=np.float64) - cx
    y = np.arange(dim, dtype=np.float64) - cy
    xx, yy = np.meshgrid(x, y)
    
    rho = np.sqrt(xx**2 + yy**2)
    theta = np.arctan2(yy, xx)
    
    if (strict is True):
        msk = (rho < rad)
    else:
        msk = (rho <= rad)

    xx[np.logical_not(msk)] = outside
    yy[np.logical_not(msk)] = outside
    rho[np.logical_not(msk)] = outside
    theta[np.logical_not(msk)] = outside

    if (normalized is True):
        xx  = xx / rad
        yy  = yy / rad
        rho = rho / rad

    if (polar is True):
        return rho, theta
    else:
        return xx, yy
    

def disc_obstructed(dim, size, obs, **kwargs):
    '''
    Create a numerical array containing a disc with central obstruction
    
    Parameters
    ----------
    dim : int
        Size of the output array
    
    size : int
        Size of the disk, representing either the radius (default) or the diameter
    
    obs : float
        Fractional size occupied by the central obstruction
    
    All other parameters are the same as for the disc() function and are described there.
    
    Returns
    -------
    disc : array
        An array containing a disc with central obstruction
    '''

    if (obs < 0) or (obs > 1):
        raise ValueError('obs value must be within [0,1]')
        return None
        
    ap_out = disc(dim, size, **kwargs)
    ap_in = disc(dim, size*obs, **kwargs)

    return ap_out - ap_in


def disc(dim, size, diameter=False, strict=False, center=(), cpix=False, invert=False, mask=False):
    '''
    Create a numerical array containing a disc.
    
    Parameters
    ----------
    dim : int
        Size of the output array
    
    size : int
        Size of the disk, representing either the radius (default) or the diameter
    
    diameter : bool, optional
        Specify if input size is the diameter. Default is 'False'
    
    strict : bool optional
        If set to Trye, size must be strictly less than (<), instead of less
        or equal (<=). Default is 'False'
    
    center : sequence, optional
        Specify the center of the disc. Default is '()', i.e. the center of the array
    
    cpix : bool optional
        If set to True, the disc is centered on pixel at position (dim//2,dim//2).
        Default is 'False', i.e. the disc is centered between 4 pixels
    
    invert : bool, optinal
        Specify if the disc must be inverted. Default is 'False'
    
    mask : bool, optional
        Specify if return value is a masked array or a numerical array. Default
        is 'False'
        
    Returns
    -------
    disc : array
        An array containing a disc with the specified parameters
    '''
    
    if (diameter is True):
        rad = size/2
    else:
        rad = size

    if (len(center) == 0):
        if (cpix is False):
            cx = (dim-1) / 2
            cy = (dim-1) / 2
        else:
            cx = dim // 2
            cy = dim // 2
    elif (len(center) == 2):
        cx = center[0]
        cy = center[1]
    else:
        raise ValueError('Error, you must pass 2 values for center')
        return None

    x = np.arange(dim, dtype=np.float64) - cx
    y = np.arange(dim, dtype=np.float64) - cy
    xx, yy = np.meshgrid(x, y)
    
    rr = np.sqrt(xx**2 + yy**2)

    if (strict is True):
        msk = (rr < rad)
    else:
        msk = (rr <= rad)
    
    if (invert is True):
        msk = np.logical_not(msk)

    if (mask is True):
        return msk
    else:
        d = np.zeros((dim, dim))
        d[msk] = 1
    
    return d


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dim  = 350
    diam = 81
    obs  = 0.2
    
    d1 = disc_obstructed(dim, diam, obs, cpix=False, center=(), strict=True, invert=True)
    
    r, t = coordinates(dim, diam, cpix=False, strict=False, center=(100, 111))

    plt.imshow(d1)
