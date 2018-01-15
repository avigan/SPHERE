# -*- coding: utf-8 -*-

'''
Aperture/pupil utility functions for pyZELDA
'''

import numpy as np
import collections

import scipy.ndimage as ndimage


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


def annulus(dim, inner_size, outer_size, **kwargs):
    '''
    Create a numerical array containing a disc with central obstruction
    
    Parameters
    ----------
    dim : int
        Size of the output array
    
    inner_size : int
        Inner size of the annulus, representing either the radius (default) or the diameter

    outer_size : int
        Outer size of the annulus, representing either the radius (default) or the diameter

    All other parameters are the same as for the disc() function and are described there.
    
    Returns
    -------
    disc : array
        An array containing a disc with central obstruction
    '''

    if (inner_size >= outer_size):
        raise ValueError('Inner radius must be smaller than outer radius')
        return None

    # outer ring
    ap_out = disc(dim, outer_size, **kwargs)

    # inner ring
    if inner_size <= 0:
        ap_final = ap_out
    else:
        ap_in = disc(dim, inner_size, **kwargs)
        ap_final = ap_out - ap_in
        
    return ap_final


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


def _rotate_interp(array, alpha, center, mode='constant', cval=0):
    '''
    Rotation around a provided center

    This is the only way to be sure where exactly is the center of rotation.

    '''
    dtype = array.dtype
    dims  = array.shape
    alpha_rad = -np.deg2rad(alpha)

    x, y = np.meshgrid(np.arange(dims[1], dtype=dtype), np.arange(dims[0], dtype=dtype))

    xp = (x-center[0])*np.cos(alpha_rad) + (y-center[1])*np.sin(alpha_rad) + center[0]
    yp = -(x-center[0])*np.sin(alpha_rad) + (y-center[1])*np.cos(alpha_rad) + center[1]

    rotated = ndimage.map_coordinates(array, [yp, xp], mode=mode, cval=cval, order=3)
    
    return rotated


def sphere_pupil(dim, diameter, dead_actuator_diameter=0.025, spiders=False, spiders_orientation=0):
    '''SPHERE pupil with dead actuators mask and spiders

    Parameters
    ----------
    dim : int
        Size of the output array
    
    diameter : int
        Diameter the disk

    dead_actuator_diameter : float
        Size of the dead actuators mask, in fraction of the pupil diameter

    spiders : bool
        Draw spiders. Default is False

    spiders_orientation : float
        Orientation of the spiders. The zero-orientation corresponds
        to the orientation of the spiders when observing in ELEV
        mode. Default is 0

    Returns
    -------
    pup : array
        An array containing a disc with the specified parameters

    '''

    # central obscuration & spiders (in fraction of the pupil)
    obs  = 0.14
    spdr = 0.005

    # main pupil
    pup = disc_obstructed(dim, diameter, obs, diameter=True, strict=False, cpix=True)

    # spiders
    spdr = max(1, spdr*dim)
    ref = np.zeros((dim, dim))
    ref[int(dim//2):, int((dim-spdr)/2+1):int((dim+spdr)//2)] = 1

    cc = dim/2
    spdr1 = _rotate_interp(ref, -5.5, (cc, cc+diameter/2))
    spdr2 = _rotate_interp(np.rot90(ref, k=1),  5.5, (cc+diameter/2, cc))
    
    spdr0 = spdr1 + spdr2
    spdr0 = np.roll(np.roll(spdr0, -1, axis=0), -1, axis=1) + np.rot90(spdr0, k=2)
    spdr0 = _rotate_interp(spdr0, 45+spiders_orientation, (cc, cc))
    spdr0 = np.roll(np.roll(spdr0, 1, axis=0), 1, axis=1)
    
    pup *= 1-spdr0
    
    # dead actuators    
    xarr = np.array([ 0.1534,  -0.0984, -0.1963,  0.2766,  0.3297])
    yarr = np.array([-0.0768,  -0.1240, -0.3542, -0.2799, -0.2799])
    for i in range(len(xarr)):
        cx = xarr[i] * diameter + diameter/2
        cy = yarr[i] * diameter + diameter/2
        
        dead = disc(dim, dead_actuator_diameter*diameter, center=(cx, cy), invert=True)

        pup *= dead

    return pup


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dim  = 350
    diam = 81
    obs  = 0.2
    
    d1 = disc_obstructed(dim, diam, obs, cpix=False, center=(), strict=True, invert=True)
    
    r, t = coordinates(dim, diam, cpix=False, strict=False, center=(100, 111))

    plt.clf()
    plt.imshow(d1)

    # plt.clf()
    # plt.imshow(sphere_pupil_internal(384, 384))
