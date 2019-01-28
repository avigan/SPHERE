# -*- coding: utf-8 -*-
'''
Images utility library

@author: avigan
'''

import collections
import numpy as np
import scipy.fftpack as fft
import scipy.ndimage as ndimage
import warnings

from astropy.convolution import convolve, Box2DKernel

####################################################################################
#
# SHIFT
#
####################################################################################


def _shift_fft(array, shift_value):
    Ndim  = array.ndim
    dims  = array.shape
    dtype = array.dtype.kind
    
    if (dtype != 'f'):
        raise ValueError('Array must be float')
    
    shifted = array
    if (Ndim == 1):
        Nx = dims[0]
        
        x_ramp = np.arange(Nx, dtype=array.dtype) - Nx//2
        
        tilt = (2*np.pi/Nx) * (shift_value[0]*x_ramp)
        
        cplx_tilt = np.cos(tilt) + 1j*np.sin(tilt)
        cplx_tilt = fft.fftshift(cplx_tilt)
        narray    = fft.fft(fft.ifft(array) * cplx_tilt)
        shifted   = narray.real
    elif (Ndim == 2):
        Nx = dims[0]
        Ny = dims[1]
        
        x_ramp = np.outer(np.full(Nx, 1.), np.arange(Ny, dtype=array.dtype)) - Nx//2
        y_ramp = np.outer(np.arange(Nx, dtype=array.dtype), np.full(Ny, 1.)) - Ny//2
        
        tilt = (2*np.pi/Nx) * (shift_value[0]*x_ramp+shift_value[1]*y_ramp)
        
        cplx_tilt = np.cos(tilt) + 1j*np.sin(tilt)        
        cplx_tilt = fft.fftshift(cplx_tilt)
        
        narray    = fft.fft2(fft.ifft2(array) * cplx_tilt)
        shifted   = narray.real
    else:
        raise ValueError('This function can shift only 1D or 2D arrays')
    
    return shifted


def _shift_interp(array, shift_value, mode='constant', cval=0):
    # Manual alternative to built-in function: slightly slower
    Ndim  = array.ndim
    dims  = array.shape
    dtype = array.dtype.kind

    if (Ndim == 1):
        pass
    elif (Ndim == 2):
        x, y = np.meshgrid(np.arange(dims[1], dtype=dtype), np.arange(dims[0], dtype=dtype))

        x -= shift_value[0]
        y -= shift_value[1]

        shifted = ndimage.map_coordinates(img, [y, x], mode=mode, cval=cval)

    return shifted


def _shift_interp_builtin(array, shift_value, mode='constant', cval=0):
    shifted = ndimage.shift(array, np.flipud(shift_value), order=3, mode=mode, cval=cval)

    return shifted


def _shift_roll(array, shift_value):
    Ndim  = array.ndim

    if (Ndim == 1):
        shifted = np.roll(array, shift_value[0])
    elif (Ndim == 2):
        shifted = np.roll(np.roll(array, shift_value[0], axis=1), shift_value[1], axis=0)
    else:
        raise ValueError('This function can shift only 1D or 2D arrays')
        
    return shifted


def shift(array, shift_value, method='fft', mode='constant', cval=0):
    '''
    Shift a 1D or 2D input array.
    
    The array can be shift either using FFT or using interpolation.

    Note that if the shifting value is an integer, the function uses
    numpy roll procedure to shift the array. The user can force to use
    np.roll, and in that case the function will round the shift value
    to the nearest integer values.
    
    Parameters
    ----------
    array : array
        The array to be shifted
    
    shift_value : float or sequence
        The shift along the axes. If a float, shift_value is the same for each axis. 
        If a sequence, shift_value should contain one value for each axis.
    
    method : str, optional
        Method for shifting the array, ('fft', 'interp'). Default is 'fft'
    
    mode : str
        Points outside the boundaries of the input are filled according 
        to the given mode ('constant', 'nearest', 'reflect' or 'wrap').
        This value is ignored if method='fft' or method='roll'. Default
        is 'constant'.
    
    cval : float, optional
        Value used for points outside the boundaries of the input if 
        mode='constant'. Default is 0
        
    Returns
    -------
    shift : array
        The shifted array

    '''

    method = method.lower()
    
    # array dimensions
    Ndim = array.ndim
    dims = array.shape
    if (Ndim != 1) and (Ndim != 2):
        raise ValueError('This function can shift only 1D or 2D arrays')

    # check that shift value is fine
    if isinstance(shift_value, collections.Iterable):
        shift_value = np.array(shift_value).ravel()
        if (shift_value.size != Ndim):
            raise ValueError('Number of dimensions in array and shift don\'t match')
    elif isinstance(shift_value, (int, float)):
        shift_value = np.full(Ndim, shift_value)
    else:
        raise ValueError('Shift value of type \'{0}\' is not allowed'.format(type(shift).__name__))    

    # check if shift values are int and automatically change method in case they are
    if (shift_value.dtype.kind == 'i'):
        method = 'roll'
    else:
        # force integer values
        if method is 'roll':
            shift_value = np.round(shift_value)
        
    # FFT limitations
    if method == 'fft':
        if np.mod(np.array(dims), 2).sum() != 0:
            raise ValueError('FFT shift only supports square images of even width')

    # detects NaN and replace them with real values
    mask = None
    nan_mask = np.isnan(array)
    if np.any(nan_mask):
        medval = np.nanmedian(array)
        array[nan_mask] = medval

        mask = np.zeros_like(array)
        mask[nan_mask] = 1
        
        mask = _shift_interp_builtin(mask, shift_value, mode='constant', cval=1)

    # shift with appropriate function                
    if (method == 'fft'):
        shifted = _shift_fft(array, shift_value)
    elif (method == 'interp'):
        shifted = _shift_interp_builtin(array, shift_value, mode=mode, cval=cval)
    elif (method == 'roll'):
        shift_value = np.round(shift_value).astype(int)
        shifted = _shift_roll(array, shift_value)
    else:
        raise ValueError('Unknown shift method \'{0}\''.format(method))

    # puts back NaN
    if mask is not None:
        shifted[mask >= 0.5] = np.nan
    
    return shifted


####################################################################################
#
# ROTATE
#
####################################################################################


def _rotate_fft(array, alpha, pad=4, x1=0, x2=0, y1=0, y2=0, cval=0):
    '''
    3 FFT shear based rotation, following Larkin et al. (1997)

    The center of rotation is exactly (dimx/2, dimy/2), i.e. at the
    center of a pixel if dimensions are even or in-between pixels if
    the dimensions are odd.
    
    Parameters
    ----------
    array : array
        The numpy array which has to be rotated
    
    alpha : float
        The rotation angle in degrees
    
    pad : int, optional 
        Padding factor. Default is 4
    
    x1,x2 : int, optional
        Borders of the original image in x
    
    y1,y2 : int, optional
        Borders of the original image in y
    
    Returns
    -------
    Rotated array

    '''

    #################################################
    # Check alpha validity and correcting if needed
    #################################################
    alpha = 1.*alpha - 360*np.floor(alpha/360)

    # FFT rotation only work in the -45:+45 range
    if alpha > 45 and alpha < 135:
        array  = np.rot90(array, k=1)
        alpha_rad = -np.deg2rad(alpha-90)
    elif alpha > 135 and alpha < 225:
        array  = np.rot90(array, k=2)
        alpha_rad = -np.deg2rad(alpha-180)
    elif alpha > 225 and alpha < 315:
        array  = np.rot90(array, k=3)
        alpha_rad = -np.deg2rad(alpha-270)
    else:
        alpha_rad = -np.deg2rad(alpha)

    ###################################
    # Preparing the frame for rotation
    ###################################
    if x1 > 0:
        px1 = (pad*array.shape[0]/2)-x1
        px2 = (pad*array.shape[0]/2)+x1
        py1 = (pad*array.shape[1]/2)-y1
        py2 = (pad*array.shape[1]/2)+y1

    pad_value = (np.array(array.shape)*(pad-1)/2).astype(np.int)
    pad_frame = np.lib.pad(array, ((pad_value[0]), (pad_value[1])), mode='constant', constant_values=cval)
    
    pad_mask = np.isnan(array)
    pad_mask = np.lib.pad(pad_mask, ((pad_value[0]), (pad_value[1])), mode='constant', constant_values=cval)
            
    # Rotate the mask, to know what part is actually the image
    pad_mask = ndimage.interpolation.rotate(pad_mask, np.rad2deg(-alpha_rad), reshape=False, order=0, mode='constant', cval=True, prefilter=False)
    
    # remove NaN from original frame
    array = np.nan_to_num(array)    
    
    # Replace part outside the image which are NaN by 0, and go into Fourier space.
    pad_frame = np.where(np.isnan(pad_frame), 0, pad_frame)
    if x1 == 0:
        pad_frame[
            np.int(((pad-1)/2)*array.shape[0]-1),
            np.int(((pad-1)/2)*array.shape[1]):np.int(((pad+1)/2)*array.shape[1])] = array[0, :] / 2
        pad_frame[
            np.int(((pad+1)/2)*array.shape[0]),
            np.int(((pad-1)/2)*array.shape[1]):np.int(((pad+1)/2)*array.shape[1])] = array[-1, :] / 2
        pad_frame[
            np.int(((pad-1)/2)*array.shape[0]):np.int(((pad+1)/2)*array.shape[0]),
            np.int(((pad-1)/2)*array.shape[1]-1)] = array[:, 0] / 2
        pad_frame[
            np.int(((pad-1)/2)*array.shape[0]):np.int(((pad+1)/2)*array.shape[0]),
            np.int(((pad+1)/2)*array.shape[1])] = array[:, -1] / 2
    elif x1 > 0:
        pad_frame[px1, py1:py2]   = pad_frame[px1, py1:py2] / 2
        pad_frame[px2-1, py1:py2] = pad_frame[px2-1, py1:py2] / 2
        pad_frame[px1:px2, py1]   = pad_frame[px1:px2, py1] / 2
        pad_frame[px1:px2, py2-1] = pad_frame[px1:px2, py2-1] / 2

    ###############################
    # Rotation in Fourier space
    ###############################
    a = np.tan(alpha_rad/2.)
    b = -np.sin(alpha_rad)

    M = -2j*np.pi*np.ones(pad_frame.shape)
    N = fft.fftfreq(pad_frame.shape[0])

    X = np.arange(-pad_frame.shape[0]/2., pad_frame.shape[0]/2.)

    pad_x   = fft.ifft((fft.fft(pad_frame, axis=0, overwrite_x=True).T * np.exp(a*((M*N).T*X).T)).T, axis=0, overwrite_x=True)
    pad_xy  = fft.ifft(fft.fft(pad_x, axis=1, overwrite_x=True) * np.exp(b*(M*X).T*N), axis=1, overwrite_x=True)
    pad_xyx = fft.ifft((fft.fft(pad_xy, axis=0, overwrite_x=True).T * np.exp(a*((M*N).T*X).T)).T, axis=0, overwrite_x=True)

    # Go back to real space
    # Put back to NaN pixels outside the image.
    pad_xyx[pad_mask] = np.NaN

    # final rotated frame
    rotated = np.real(pad_xyx[
        np.int(((pad-1)/2)*array.shape[0]):np.int(((pad+1)/2)*array.shape[0]),
        np.int(((pad-1)/2)*array.shape[1]):np.int(((pad+1)/2)*array.shape[1])]).copy()
    
    return rotated


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


def _rotate_interp_builtin(array, alpha, center, mode='constant', cval=0):
    '''
    Rotation with the built-in rotation function.

    The center of rotation is exactly ((dimx-1) / 2, (dimy-1) / 2),
    whatever the odd/even-ity of the image dimensions

    '''
    alpha_rad = -np.deg2rad(alpha)

    # center of rotation: ((dims[1]-1) / 2, (dims[0]-1) / 2)
    rotated = ndimage.interpolation.rotate(array, np.rad2deg(-alpha_rad), reshape=False, order=3, mode=mode, cval=cval)

    return rotated


def _rotate_roll(array, alpha):
    '''
    Roll rotation for multiples of 90 deg

    The center of rotation is exactly ((dimx-1) / 2, (dimy-1) / 2),
    whatever the odd/even-ity of the image dimensions

    '''
    k = (alpha / 90) % 4
    
    rotated = np.rot90(array, k=k)
    
    return rotated


def rotate(array, value, center=None, method='interp', mode='constant', cval=0):
    '''
    Clockwise rotation of images

    Parameters
    ----------
    array : array
        The numpy array which has to be rotated
    
    value : float
        The rotation angle in degrees
    
    center : array
        Center of rotation. Default value is (dims[1] // 2,
        dims[0] // 2), i.e. always at the center of a pixel.
    
    method : str, optional    
        Method for shifting the array, ('fft', 'interp'). Default is
        'interp'

    mode : str    
        Points outside the boundaries of the input are filled
        according to the given mode ('constant', 'nearest', 'reflect'
        or 'wrap').  This value is ignored if method='fft'. Default is
        'constant'.
    
    cval : float, optional
        Value used for points outside the boundaries of the input if 
        mode='constant'. Default is 0.0
        
    Returns
    -------
    rotated : array
        The rotated array

    '''

    method = method.lower()

    array = array.copy()
    
    # array dimensions
    Ndim = array.ndim
    dims = array.shape
    if (Ndim != 2):
        raise ValueError('This function can rotate only 2D arrays')

    # center of rotation
    if center is None:
        center = np.array((dims[1] // 2, dims[0] // 2))
    else:
        center = np.array(center).ravel()
        if (center.size != 2):
            raise ValueError('You must pass 2 values for center')
        
        if method == 'fft':
            warnings.warn('\'center\' is ignored with method=\'fft\'. ' +
                          'The center of rotation is exactly at (dimx/2, dimy/2).')

    # check rotation value
    if not isinstance(value, (int, float)):
        raise ValueError('Rotation value of type \'{0}\' is not allowed'.format(type(value).__name__))
            
    # clockwise
    value *= -1

    # check if shift values are int and automatically change method in case they are
    if ((value % 90) == 0):
        method = 'roll'

    # FFT limitations
    if method == 'fft':
        if np.mod(np.array(dims), 2).sum() != 0:
            raise ValueError('FFT rotate only supports square images of even width')
    
    # detects NaN and replace them with real values
    mask = None
    nan_mask = np.isnan(array)
    if np.any(nan_mask):
        medval = np.nanmedian(array)
        array[nan_mask] = medval

        mask = np.zeros_like(array)
        mask[nan_mask] = 1
        
        mask = _rotate_interp(mask, value, center, mode='constant', cval=1)

    # rotate with appropriate function
    if (method == 'fft'):
        rotated = _rotate_fft(array, value)
    elif (method == 'interp'):
        rotated = _rotate_interp(array, value, center, mode=mode, cval=cval)
    elif (method == 'roll'):
        rotated = _rotate_roll(array, value)
    else:
        raise ValueError('Unknown rotate method \'{0}\''.format(method))

    # puts back NaN
    if mask is not None:
        rotated[mask >= 0.5] = np.nan
    
    return rotated


####################################################################################
#
# SCALE
#
####################################################################################


def _fft_floating_origin(array, center_direct=None, center_fourier=None, inverse=False, cc=False, ee=False, ce=False, ec=False):
    '''
    FFT calculation with floating origin

    Parameters
    ----------
    array : array_like
        The array to be transformed

    center_direct : array_like
        Center in direct space

    center_fourier : array_like
        Center in Fourier space
    '''

    dim   = array.shape[0]
    dtype = array.dtype

    # direction of the FFT
    if inverse:
        sign = 1
    else:
        sign = -1
    
    # centers
    if center_direct is None:
        center_direct = np.array(((dim-1)/2, (dim-1)/2))

    if center_fourier is None:
        center_fourier = np.array(((dim-1)/2, (dim-1)/2))

    if cc:
        center_direct  = np.array(((dim)/2, (dim)/2))
        center_fourier = np.array(((dim)/2, (dim)/2))

    if ce:
        center_direct  = np.array(((dim)/2, (dim)/2))
        center_fourier = np.array(((dim-1)/2, (dim-1)/2))

    if ec:
        center_direct  = np.array(((dim-1)/2, (dim-1)/2))
        center_fourier = np.array(((dim)/2, (dim)/2))

    if ee:
        center_direct  = np.array(((dim-1)/2, (dim-1)/2))
        center_fourier = np.array(((dim-1)/2, (dim-1)/2))

    # usefull parameters
    x, y = np.meshgrid(np.arange(dim, dtype=dtype), np.arange(dim, dtype=dtype))

    # shift in Fourier space by multiplying in direct space
    tmp = array * np.exp((-sign) * 2*np.pi*1j * (center_fourier[0]*x + center_fourier[1]*y) / dim)
    if sign < 0:
        array_f = fft.fft2(tmp) / dim**2
    else:
        array_f = fft.ifft2(tmp) * dim**2

    # shift in direct space by multiplying in Fourier space
    array_f *= np.exp((-sign) * 2*np.pi*1j * (center_direct[0]*x + center_direct[1]*y) / dim)

    # normalisation
    array_f *= np.exp(sign * 2*np.pi*1j / dim * (center_direct*center_fourier).sum())

    return array_f

    
def _scale_fft(array, scale_value, alt_criterion=False):
    dim   = array.shape[0]    # square-even images checked before
    dtype = array.dtype.kind

    # take the first of the two values. Check done before that the two
    # values are identical
    zoom_io = scale_value[0]

    # We want Fechs *close* *to* N'/N" where N" = N + 2*KF, N' = N + 2*KD
    #   => N" = 2*round(N'/(2*Fechs))
    #   => KF = (N"-N)/2 = round(N'/(2*Fechs) - N/2) = round(N/2*(1/Fechs-1) + KD/2)
    # We call yy=N/2*(1/Fechs-1) +KD/2
    # We minimize this difference between the `ideal' N" and its closest integer value      
    # Compared to the ALTernate criterion below, this one favors small
    # values of N" i.e. little truncation in Fourier space.  
    kd_array = np.arange(dim/2 + 1, dtype=np.int)
    yy = dim/2 * (zoom_io - 1) + kd_array.astype(np.float)*zoom_io
    kf_array = np.round(yy).astype(np.int)

    tmp = np.abs(yy-kf_array)
    # tmp[tmp == 0] = np.nan
    imin = np.nanargmin(tmp)

    kd_io = kd_array[imin]
    kf_io = kf_array[imin]

    # alternate criterion: minimize error |Fechs-N'/N"|
    if alt_criterion:
        error = np.abs(1/zoom_io - (dim + 2*kd_array) / (dim + 2*kf_array))
        iminerror = np.argmin(error)

        kd_io = kd_array[iminerror]
        kf_io = kf_array[iminerror]
    
    # Extract a part of, or expand, array to dim_p pixels
    dim_p = int(dim + 2*kd_io)
    
    if dim_p > dim:
        # kd is positive
        tmp = np.zeros((dim_p, dim_p), dtype=dtype)
        tmp[kd_io:kd_io+dim, kd_io:kd_io+dim] = array
    else:
        # kd is negative
        tmp = array[-kd_io:-kd_io+dim_p, -kd_io:-kd_io+dim_p]

    # Fourier-transform the result
    array_f = _fft_floating_origin(tmp, ec=True) * dim_p**2
    
    # Extract a part of or expand the result to dim_pp pixels
    dim_pp = int(dim + 2*kf_io)
    
    if dim_pp > dim_p:
        tmp = np.zeros((dim_pp, dim_pp), dtype=np.complex)
        tmp[(dim_pp-dim_p)//2:(dim_pp+dim_p)//2, (dim_pp-dim_p)//2:(dim_pp+dim_p)//2] = array_f
    else:
        tmp = array_f[kd_io-kf_io:kd_io-kf_io+dim_pp, kd_io-kf_io:kd_io-kf_io+dim_pp]

    # inverse Fourier-transform the result
    tmp = _fft_floating_origin(tmp, inverse=True, ce=True)
    array_scaled = tmp.real  / dim_pp**2
    del tmp

    # Extract a part of or expand the result back to NP pixels
    if dim_pp > dim:
        # kf is positive
        scaled = array_scaled[kf_io:kf_io+dim, kf_io:kf_io+dim]
    else:
        # kf is negative
        scaled = array*0
        scaled[-kf_io:-kf_io+dim_pp, -kf_io:-kf_io+dim_pp] = array_scaled

    return scaled
    

def _scale_interp(array, scale_value, center, mode='constant', cval=0):
    Ndim  = array.ndim
    dims  = array.shape
    dtype = array.dtype.kind

    if (Ndim == 1):
        pass
    elif (Ndim == 2):
        x, y = np.meshgrid(np.arange(dims[1], dtype=dtype), np.arange(dims[0], dtype=dtype))

        nx = (x - center[0]) / scale_value[0] + center[0]
        ny = (y - center[1]) / scale_value[1] + center[1]
        
        scaled = ndimage.map_coordinates(array, [ny, nx], mode=mode, cval=cval)

    return scaled


def _scale_interp_builtin(array, scale_value, mode='constant', cval=0):
    scaled = ndimage.zoom(array, scale_value, order=3, mode=mode, cval=cval)

    return scaled
    

def scale(array, scale_value, center=None, new_dim=None, method='fft', mode='constant', cval=0):
    '''
    Rescale a 2D input array.
    
    The array can be rescaled either using FFT or interpolation. When
    using interpolation, the scaling factor can be different in x and
    y to enable anamorphism correction.
    
    Parameters
    ----------
    array : array
        The array to be shifted
    
    scale_value : float or sequence
        The scaling along the axes. If a float, scale_value is the same for each axis.
        If a sequence, scale_value should contain one value for each axis.

    center : array
        Center of the scaling. Default value is (dims[1] // 2, dims[0] // 2),
        i.e. always at the center of a pixel.

    new_dim : array
        Directly specifies the new dimensions of the resulting array. If provided,
        the 'scale_value' and 'center' parameters are ignored.
    
    method : str, optional
        Method for shifting the array, ('fft', 'interp'). Default is 'fft'
    
    mode : str    
        Points outside the boundaries of the input are filled
        according to the given mode ('constant', 'nearest', 'reflect'
        or 'wrap').  This value is ignored if method='fft'. Default is
        'constant'.
    
    cval : float, optional
        Value used for points outside the boundaries of the input if 
        mode='constant'. Default is 0

    Returns
    -------
    shift : array
        The shifted array

    '''

    method = method.lower()

    # no scaling needed
    if scale_value == 1:
        return array
    
    # array dimensions
    Ndim  = array.ndim
    dims = array.shape
    if (Ndim != 2):
        raise ValueError('This function can scale only 2D arrays')

    # center of scaling
    if center is None:
        center = np.array((dims[1] // 2, dims[0] // 2))
    else:
        center = np.array(center).ravel()
        if (center.size != 2):
            raise ValueError('You must pass 2 values for center')

        if method == 'fft':
            warnings.warn('Center is not taken into account with method=\'fft\'. ' +
                          'The center of rotation is exactly at ((dimx-1)/2, dimy/2).')
    
    # check that scale value is fine
    if isinstance(scale_value, collections.Iterable):
        scale_value = np.array(scale_value).ravel()
        if (scale_value.size != Ndim):
            raise ValueError('Number of dimensions in array and scale value don\'t match')
    elif isinstance(scale_value, (int, float)):
        scale_value = np.full(Ndim, scale_value)
    else:
        raise ValueError('Shift value of type \'{0}\' is not allowed'.format(type(scale_value).__name__))

    # user-specified dimensions
    if new_dim is not None:
        new_dim = np.array(new_dim).ravel()
        if (new_dim.size != 2):
            raise ValueError('You must pass 2 values for new_dim')

        warnings.warn('\'center\' and \'scale_value\' are ignored when \'new_dim\' is specified.')

        if method == 'fft':
            warnings.warn('method = \'fft\' is not supported when \'new_dim\' is specified.' +
                          'Switching to method = \'interp\'.')

        method = 'interp_builtin'
        scale_value = np.flipud(new_dim) / np.array(dims)

    # FFT limitations
    if method == 'fft':
        if (np.mod(np.array(dims), 2).sum() != 0) or (dims[0] != dims[1]):
            raise ValueError('FFT scale only supports square images of even width')

        if scale_value[0] != scale_value[1]:
            raise ValueError('FFT scale only supports identical factor along the two dimensions')
        
    # detects NaN and replace them with real values
    mask = None
    nan_mask = np.isnan(array)
    if np.any(nan_mask):
        medval = np.nanmedian(array)
        array[nan_mask] = medval

        mask = np.zeros_like(array)
        mask[nan_mask] = 1

        if (method == 'fft'):
            mask = _scale_interp(mask, scale_value, center, mode='constant', cval=1)            
        elif (method == 'interp'):
            mask = _scale_interp(mask, scale_value, center, mode='constant', cval=1)
        elif (method == 'interp_builtin'):
            mask = _scale_interp_builtin(array, scale_value, mode='constant', cval=1)

    # scale with appropriate function
    if (method == 'fft'):
        scaled = _scale_fft(array, scale_value)
    elif (method == 'interp'):
        scaled = _scale_interp(array, scale_value, center, mode=mode, cval=cval)
    elif (method == 'interp_builtin'):
        scaled = _scale_interp_builtin(array, scale_value, mode=mode, cval=cval)
    else:
        raise ValueError('Unknown scale method \'{0}\''.format(method))

    # puts back NaN
    if mask is not None:
        scaled[mask >= 0.5] = np.nan
    
    return scaled
    

####################################################################################
#
# CLEANING
#
####################################################################################


def sigma_filter(img, box=5, nsigma=3, iterate=False, return_mask=False, max_iter=20, _iters=0, _mask=None):
    '''
    Performs sigma-clipping over an image

    Adapted from the IDL function with the same name in the astron library.
    
    Parameters
    ----------
    img : array
        The input image
    
    box : int, optional
        Box size for the sigma-clipping. Default is 5 pixel
    
    nsigma : float, optional
        Sigma value. Default if 3.
    
    iterate : bool, optional
        Controls if the filtering is iterative. Default is False

    return_mask : bool
        If True, returns a mask to identify the clipped values. Default is False
    
    max_iter : int, optional
        Maximum number of iterations. Default is 20
    
    _iters : int (internal)
        Internal counter to keep track during iterative sigma-clipping

    _mask : array_like (internal)
        Keep track of bad pixels over the iterations
    
    Returns
    -------
    return_value : array
        Input image with clipped values
    
    '''

    # clip bad pixels
    box2 = box**2

    kernel = Box2DKernel(box)
    img_clip = (convolve(img, kernel)*box2 - img) / (box2-1)

    imdev = (img - img_clip)**2
    fact = nsigma**2 / (box2-2)
    imvar = fact*(convolve(imdev, kernel)*box2 - imdev)

    # following solution is faster but does not support bad pixels
    # see avigan/VLTPF#49
    # img_clip = (ndimage.uniform_filter(img, box, mode='constant')*box2 - img) / (box2-1)

    # imdev = (img - img_clip)**2
    # fact = nsigma**2 / (box2-2)
    # imvar = fact*(ndimage.uniform_filter(imdev, box, mode='constant')*box2 - imdev)
    
    wok = np.nonzero(imdev < imvar)
    nok = wok[0].size
    
    # copy good pixels in clipped image
    if (nok > 0):
        img_clip[wok] = img[wok]

    # create _mask at first iteration
    if _mask is None:
        _mask = np.zeros_like(img, dtype=np.bool)

    # identify clipped pixels
    _mask[img != img_clip] = True

    # iterations
    nchange = img.size - nok
    if (iterate is True):
        _iters = _iters+1
        if (_iters >= max_iter) or (nchange == 0):
            if return_mask:
                return img_clip, _mask
            else:
                return img_clip

        return sigma_filter(img_clip, box=box, nsigma=nsigma, iterate=iterate,
                            return_mask=return_mask, _iters=_iters, _mask=_mask)
        
    if return_mask:
        return img_clip, _mask
    else:
        return img_clip


def fix_badpix_vip(img, bpm, box=5):
    '''
    Corrects the bad pixels, marked in the bad pixel mask.

    The bad pixels are replaced by the median of the adjacent pixels
    in a box of the povided box size. This function is very fast but
    works best with isolated (sparse) pixels or very small clusters.

    Copied and adapted from the the Vortex Image Processing package,
    https://github.com/vortex-exoplanet/VIP, in which the function is
    called fix_badpix_isolated.

    This version is improved with respect to the VIP one by replacing
    the bad pixels with NaNs in the image before applying the
    median_filter. This allows to make sure that adjacent bad pixels
    will not be taken into account when calculating the median.
    
    Parameters
    ----------
    img : array_like
        Input 2D image
    
    bpm : array_like, optional
        Input bad pixel map. Good pixels have a value of 0, bad pixels
        a value of 1.

    box : odd int, optional
        The size the box (box x box) of adjacent pixels for the
        median filter. Default value is 5
    
    Return
    ------
    img_clean : array_like
        Cleaned image

    '''
    
    if not img.ndim == 2:
        raise ValueError('Main input is not a 2D array')
    
    if not bpm.ndim == 2:
        raise ValueError('Bad pixel map input is not a 2D array')
    
    if box % 2 == 0:
        raise ValueError('Box size of the median blur kernel must be an odd integer')

    bpm = bpm.astype('bool')

    bp = np.where(bpm)
    
    img_clean = img.copy()
    img_clean[bp] = np.nan

    smoothed = ndimage.median_filter(img_clean, box, mode='mirror')
    img_clean[bp] = smoothed[bp]

    # replace uncorrected bad pixels with original value
    mask = ~np.isfinite(img_clean)
    img_clean[mask] = img[mask]
    
    return img_clean


def fix_badpix(img, bpm, npix=8, weight=False):
    '''Corrects the bad pixels, marked in the bad pixel mask.

    It will fill in bad pixels by finding the NPIX nearest good
    pixels, toss the highest and lowest ones of the bunch, and then
    arithmatically average. Additional it will weight adjacent pixels
    by inverse of their distances in averaging process if the option
    is selected.

    Important warning: to make computation faster, the weighing is not
    applied for bad pixels located within a few pixels from the edges
    of the image.

    Parameters
    ----------
    img : array_like
        Input 2D image
    
    bpm : array_like, optional    
        Input bad pixel map. Good pixels have a value of 0, bad pixels
        a value of 1.

    npix : int, optional    
        The number of adjacent good pixels used for the estimation bad
        pixel value. Default value is 8

    weight : bool, optional
        Weigh good pixel by inverse of their distance in the averaging
        process. Default is False
    
    Return
    ------
    img_clean : array_like
        Cleaned image

    '''
    # new arrays
    img = img.copy()
    bpm = (bpm != 0)

    # bad pixels
    bp  = np.where(bpm)
    nbp = bp[0].size    
    if nbp == 0:
        return img

    # usefull parameters
    ddmin = 2
    ddmax = 100
    shape = img.shape

    # create default distance array
    dd = ddmin
    xx, yy = np.meshgrid(np.arange(2*dd+1)-dd, np.arange(2*dd+1)-dd)
    dist_default = np.sqrt(xx**2 + yy**2)

    bpm = np.logical_not(bpm)
    for cx, cy in zip(bp[1], bp[0]):
        # default search box is 2*dd+1 pixel
        dd = ddmin

        # determine search region
        found = False
        while not found:
            x0 = max(cx-dd, 0)
            x1 = min(cx+dd+1, shape[-1])
            y0 = max(cy-dd, 0)
            y1 = min(cy+dd+1, shape[-2])

            bpm_sub = bpm[y0:y1, x0:x1]
            img_sub = img[y0:y1, x0:x1]
            
            if bpm_sub.sum() < npix:
                dd = dd + 2
            else:
                found = True

            if dd > ddmax:
                break
            
        # distance to adjacent good pixels
        if dd == ddmin:
            # get default array if dd unchanged
            dist = dist_default
        else:
            # otherwise recompute one
            xx, yy = np.meshgrid(np.arange(2*dd+1)-dd, np.arange(2*dd+1)-dd)
            dist = np.sqrt(xx**2 + yy**2)

        # no weighing if we at the edges
        if (bpm_sub.shape != (2*dd+1, 2*dd+1)):
            dist = np.ones_like(bpm_sub)

        # keep good pixels
        good_pix  = img_sub[bpm_sub]
        good_dist = dist[bpm_sub]
        
        # sort them by distance
        ii = np.argsort(good_dist)
        good_pix  = good_pix[ii]
        good_dist = good_dist[ii]

        # get values of relevant pixels
        mm = np.where(good_dist <= good_dist[npix-1])
        good_pix  = good_pix[mm]
        good_dist = good_dist[mm]

        ii = np.argsort(good_pix)
        good_pix  = good_pix[ii]
        good_dist = good_dist[ii]
        
        # calculate new pixel value, tossing the highest and lowest
        # pixels of the bunch, then weighting by the inverse of the
        # distances if desired
        if weight:
            final_dist = good_dist[1:-1]
            new_val = np.sum(good_pix[1:-1] / final_dist)
            new_val = new_val / np.sum(1/final_dist)
        else:
            new_val = np.mean(good_pix[1:-1])
            
        img[cy, cx] = new_val
            
    return img

####################################################################################
#
# PROFILES
#
####################################################################################


def profile(img, type='mean', step=1, mask=None, center=None, rmax=0, clip=True, exact=False):
    '''
    Azimuthal statistics of an image

    Parameters
    ----------
    img : array
        Image on which the profiles
        
    ptype : str, optional
        Type of profile. Allowed values are mean, std, var, median, min, max. Default is mean.
    
    mask : array, optional
        Mask for invalid values (must have the same size as image)
        
    center : array_like, optional
        Center of the image

    rmax : float
        Maximum radius for calculating the profile, in pixel. Default is 0 (no limit)
    
    clip : bool, optional
        Clip profile to area of image where there is a full set of data
        
    exact : bool, optional
        Performs an exact estimation of the profile. This can be very long for 
        large arrays. Default is False, which rounds the radial distance to the 
        closest 1 pixel.
    
    Returns
    -------
    prof : array
        1D profile vector
        
    rad : array
        Separation vector, in pixel
    '''
    
    # array dimensions
    dimx = img.shape[1]
    dimy = img.shape[0]

    # center
    if center is None:
        center = (dimx // 2, dimy // 2)

    # masking
    if mask is not None:
        # check size
        if mask.shape != img.shape:
            raise ValueError('Image and mask don''t have the same size. Returning.')

        img[mask == 0] = np.nan
        
    # intermediate cartesian arrays
    x = np.arange(dimx, dtype=np.int64) - center[0]
    y = np.arange(dimy, dtype=np.int64) - center[1]
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    
    # rounds for faster calculation
    if not exact:
        rr = np.round(rr, decimals=0)
    
    # find unique radial values
    uniq = np.unique(rr, return_inverse=True, return_counts=True)
    r_uniq_val = uniq[0]
    r_uniq_inv = uniq[1]
    r_uniq_cnt = uniq[2]

    # number of elements
    if clip:
        extr  = np.abs(np.array((x[0], x[-1], y[0], y[-1])))
        r_max = extr.min()
        i_max = int(r_uniq_val[r_uniq_val <= r_max].size)
    else:
        r_max = r_uniq_val.max()
        i_max = r_uniq_val.size

    # limit extension of profile
    if (rmax > 0):
        r_max = rmax
        i_max = int(r_uniq_val[r_uniq_val <= r_max].size)
        
    t_max = r_uniq_cnt[0:i_max].max()

    # intermediate polar array
    polar = np.empty((i_max, t_max), dtype=img.dtype)
    polar.fill(np.nan)
    
    img_flat = img.ravel()
    for r in range(i_max):
        cnt = r_uniq_cnt[r]
        val = img_flat[r_uniq_inv == r]
        polar[r, 0:cnt] = val
            
    # calculate profile
    rad  = r_uniq_val[0:i_max]

    type = type.lower()
    if step == 1:
        # fast statistics if step=1
        if type == 'mean':
            prof = np.nanmean(polar, axis=1)
        elif type == 'std':
            prof = np.nanstd(polar, axis=1, ddof=1)
        elif type == 'var':
            prof = np.nanvar(polar, axis=1)
        elif type == 'median':
            prof = np.nanmedian(polar, axis=1)
        elif type == 'min':
            prof = np.nanmin(polar, axis=1)
        elif type == 'max':
            prof = np.nanmax(polar, axis=1)
        else:
            raise ValueError('Unknown statistics ptype = {0}. Allowed values are mean, std, var, median, min and max'.format(type))
    else:
        # slower if we need step > 1
        prof = np.zeros(i_max, dtype=img.dtype)
        for r in range(i_max):
            idx = ((rad[r]-step/2) <= rad) & (rad <= (rad[r]+step/2))
            val = polar[idx, :]
            
            if type == 'mean':
                prof[r] = np.nanmean(val)
            elif type == 'std':
                prof[r] = np.nanstd(val)
            elif type == 'var':
                prof[r] = np.nanvar(val)
            elif type == 'median':
                prof[r] = np.nanmedian(val)
            elif type == 'min':
                prof[r] = np.nanmin(val)
            elif type == 'max':
                prof[r] = np.nanmax(val)
            else:
                raise ValueError('Unknown statistics ptype = {0}. Allowed values are mean, std, var, median, min and max'.format(type))

    return prof, rad


####################################################################################
#
# FILTERING
#
####################################################################################

def median(img, dim):
    '''
    Apply median filtering to an image 

    Parameters
    ----------
    img : array
        Image on which the profiles
        
    dim : int
        Size of the median filter
    
    Returns
    -------
    img_filt : array
        Median filtered image
    '''
    img_filt = img - ndimage.filters.median_filter(img, size=(dim, dim))
    
    return img_filt
    

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from vigan.utils import aperture as aper
    from astropy.io import fits

    #
    # 1D - shift
    #
    s = 12.5
    x = np.arange(300)
    
    cos = np.cos(x/299*4*np.pi)
    cos_roll   = shift(cos, shift_value=s, method='roll')
    cos_interp = shift(cos, shift_value=s, method='interp')
    cos_fft    = shift(cos, shift_value=s, method='fft')

    fig = plt.figure('Shift - 1D', figsize=(12, 10))
    plt.clf()    
    ax = fig.add_subplot(111)
    
    ax.plot(x, cos, marker='', label='Original')
    ax.plot(x, cos_roll, marker='', label='Roll', linewidth=8)
    ax.plot(x, cos_interp, marker='', label='Interp', linewidth=4)
    ax.plot(x, cos_fft, marker='', label='FFT', linewidth=1)

    ax.legend()
    
    plt.tight_layout()
    
    #
    # 2D - shift
    #
    c = (75, 200)
    s = (80.5, -60.2)

    # images
    img = aper.disc(300, 50, center=c, diameter=True)    
    img_roll   = shift(img, shift_value=s, method='roll')
    img_interp = shift(img, shift_value=s, method='interp')
    img_fft    = shift(img, shift_value=s, method='fft')

    # plot
    fig = plt.figure('Shift - 2D', figsize=(12, 12))
    plt.clf()

    ax = fig.add_subplot(221)
    ax.imshow(img, aspect='equal')
    ax.plot(c[0], c[1], marker='+')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Original')

    ax = fig.add_subplot(222)
    ax.imshow(img_roll, aspect='equal')
    ax.plot(c[0]+s[0], c[1]+s[1], marker='+')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Roll')

    ax = fig.add_subplot(223)
    ax.imshow(img_interp, aspect='equal')
    ax.plot(c[0]+s[0], c[1]+s[1], marker='+')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Interp')

    ax = fig.add_subplot(224)
    ax.imshow(img_fft, aspect='equal')
    ax.plot(c[0]+s[0], c[1]+s[1], marker='+')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('FFT')

    plt.tight_layout()

    #
    # 2D - rotate
    #
    c = (75, 200)
    v = -10.5

    # images
    img = aper.disc(300, 50, center=c, diameter=True)
    
    img_roll   = rotate(img, value=v, method='roll')
    img_interp = rotate(img, value=v, method='interp')
    img_fft    = rotate(img, value=v, method='fft')

    # plot
    fig = plt.figure('Rotate', figsize=(12, 12))
    plt.clf()

    ax = fig.add_subplot(221)
    ax.imshow(img, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Original')

    ax = fig.add_subplot(222)
    ax.imshow(img_roll, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Roll')

    ax = fig.add_subplot(223)
    ax.imshow(img_interp, aspect='equal', vmin=0, vmax=0.1)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Interp')

    ax = fig.add_subplot(224)
    ax.imshow(img_fft, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('FFT')
    
    plt.tight_layout()
        
    #
    # 2D - scale
    #
    c = (80, 220)
    s = 0.91907

    # images
    img = aper.disc(300, 50, center=c, diameter=True)
    img_interp = scale(img, scale_value=s, center=(300//2, 300/2), method='interp')
    img_interp_newdim = scale(img, 1, new_dim=(250, 75), method='interp')
    img_fft = scale(img, scale_value=s, center=(300//2, 300/2), method='fft')
    
    # plot
    fig = plt.figure('Shift - 2D', figsize=(12, 12))
    plt.clf()

    ax = fig.add_subplot(221)
    ax.imshow(img, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Original')

    ax = fig.add_subplot(222)
    ax.imshow(img_interp, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Interp')

    ax = fig.add_subplot(223)
    ax.imshow(img_interp_newdim, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('Interp - new_dim')

    ax = fig.add_subplot(224)
    ax.imshow(img_fft, aspect='equal')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title('FFT')

    plt.tight_layout()
    
