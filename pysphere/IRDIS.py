import os
import glob
import pandas as pd
import subprocess
import numpy as np
import astropy.coordinates as coord
import astropy.units as units
import scipy.ndimage as ndimage
import scipy.interpolate as interp
import scipy.optimize as optim
import shutil
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors

from astropy.io import fits
from astropy.time import Time
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages

import pysphere.utils.imutils as imutils
import pysphere.utils.aperture as aperture
import pysphere.transmission as transmission
import pysphere.ReductionPath as ReductionPath


# keywords to be saved
keywords = [
    # standard
    'INSTRUME', 
    'OBJECT', 'DATE-OBS', 'DATE', 'HIERARCH ESO DET FRAM UTC',

    # DPR
    'HIERARCH ESO DPR CATG', 'HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR TECH',
        
    # coordinates
    'HIERARCH ESO TEL GEOLAT', 'HIERARCH ESO TEL GEOLON', 'HIERARCH ESO TEL GEOELEV',
    'HIERARCH ESO INS4 DROT2 RA', 'HIERARCH ESO INS4 DROT2 DEC',
    'HIERARCH ESO TEL ALT', 'HIERARCH ESO TEL AZ',

    # SAXO
    'HIERARCH ESO AOS TTLOOP STATE', 'HIERARCH ESO AOS HOLOOP STATE',
    'HIERARCH ESO AOS IRLOOP STATE', 'HIERARCH ESO AOS PUPLOOP STATE',
    'HIERARCH ESO AOS VISWFS MODE',

    # CPI
    'HIERARCH ESO SEQ ARM',
    'HIERARCH ESO INS COMB ICOR', 'HIERARCH ESO INS COMB IFLT', 'HIERARCH ESO INS COMB POLA',
    'HIERARCH ESO INS4 FILT2 NAME',
    'HIERARCH ESO INS4 DROT2 BEGIN', 'HIERARCH ESO INS4 DROT2 END',
    'HIERARCH ESO INS4 DROT2 POSANG', 'HIERARCH ESO INS4 DROT2 MODE', 
    
    # IFS
    'HIERARCH ESO INS2 MODE', 'HIERARCH ESO INS2 COMB IFS',

    # IRDIS
    'HIERARCH ESO INS1 MODE',
    'HIERARCH ESO INS1 FILT NAME', 'HIERARCH ESO INS1 OPTI2 NAME',
    'HIERARCH ESO INS1 PAC X', 'HIERARCH ESO INS1 PAC Y',
    
    # detector
    'HIERARCH ESO DET SEQ1 DIT', 'HIERARCH ESO DET NDIT',

    # observing conditions
    'HIERARCH ESO TEL AIRM START', 'HIERARCH ESO TEL AIRM END',
    'HIERARCH ESO TEL AMBI FWHM START', 'HIERARCH ESO TEL AMBI FWHM END', 'HIERARCH ESO TEL IA FWHM',
    'HIERARCH ESO TEL AMBI TAU0', 'HIERARCH ESO TEL AMBI TEMP',
    'HIERARCH ESO TEL AMBI WINDSP', 'HIERARCH ESO TEL AMBI WINDDIR'
]

# short keywords
keywords_short = keywords.copy()
for idx in range(len(keywords_short)):
    key = keywords_short[idx]
    if key.find('HIERARCH ESO ') != -1:
        keywords_short[idx] = key[13:]
        
# useful parameters
nwave = 2
pixel = 12.25

wave_cal_lasers = [0.9877, 1.1237, 1.3094, 1.5451]

# specify for each recipe which other recipes need to have been executed before
recipe_requirements = {
    'sort_files': [],
    'sort_frames': ['sort_files'],
    'check_files_association': ['sort_files'],
    'sph_ird_cal_dark': ['sort_files'],
    'sph_ird_cal_detector_flat': ['sort_files'],
    'sph_ird_preprocess_science': ['sort_files', 'sort_frames', 'sph_ird_cal_dark', 'sph_ird_cal_detector_flat'],
    'sph_ird_star_center': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science'],
    'sph_ird_combine_data': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science', 'sph_ird_star_center']
}


def check_recipe_execution(recipe_execution, recipe_name):
    '''
    Check execution of previous recipes for a given recipe.

    Parameters
    ----------
    recipe_execution : dict
        Status of executed recipes

    recipe_name : str
        Name of the current recipe

    Returns
    -------
    execute_recipe : bool
        Current recipe can be executed safely
    '''
    requirements = recipe_requirements[recipe_name]

    execute_recipe = True
    missing = []
    for r in requirements:
        if not recipe_execution[r]:
            execute_recipe = False
            missing.append(r)

    if not execute_recipe:
        raise ValueError('{0} cannot executed because some files have been '.format(recipe_name) +
                         'removed from the reduction directory ' +
                         'or the following recipes have not been executed: {0}. '.format(missing))

    return execute_recipe

    
def parallatic_angle(ha, dec, geolat):
    '''
    Parallactic angle of a source in degrees

    Parameters
    ----------
    ha : array_like
        Hour angle, in hours

    dec : float
        Declination, in degrees

    geolat : float
        Observatory declination, in degrees

    Returns
    -------
    pa : array_like
        Parallactic angle values
    '''
    pa = -np.arctan2(-np.sin(ha),
                     np.cos(dec) * np.tan(geolat) - np.sin(dec) * np.cos(ha))

    if (dec >= geolat):
        pa[ha < 0] += 360*units.degree
    
    return np.degrees(pa)


def compute_times(frames_info):
    '''
    Compute the various timestamps associated to frames

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames
    '''

    # get necessary values
    time_start = frames_info['DATE-OBS'].values
    time_end   = frames_info['DET FRAM UTC'].values
    time_delta = (time_end - time_start) / frames_info['DET NDIT'].values.astype(np.int)
    DIT        = np.array(frames_info['DET SEQ1 DIT'].values.astype(np.float)*1000, dtype='timedelta64[ms]')

    # calculate UTC time stamps
    idx = frames_info.index.get_level_values(1).values
    ts_start = time_start + time_delta * idx
    ts       = time_start + time_delta * idx + DIT/2
    ts_end   = time_start + time_delta * idx + DIT

    # calculate mjd
    geolon = coord.Angle(frames_info['TEL GEOLON'].values[0], units.degree)
    geolat = coord.Angle(frames_info['TEL GEOLAT'].values[0], units.degree)
    geoelev = frames_info['TEL GEOELEV'].values[0]

    utc = Time(ts_start.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    mjd_start = utc.mjd
    
    utc = Time(ts.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    mjd = utc.mjd

    utc = Time(ts_end.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    mjd_end = utc.mjd
    
    # update frames_info
    frames_info['TIME START'] = ts_start
    frames_info['TIME']       = ts
    frames_info['TIME END']   = ts_end

    frames_info['MJD START']  = mjd_start
    frames_info['MJD']        = mjd
    frames_info['MJD END']    = mjd_end


def compute_angles(frames_info):
    '''
    Compute the various angles associated to frames: RA, DEC, parang,
    pupil offset, final derotation angle

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames
    '''

    # derotator drift check and correction
    date_fix = Time('2016-07-12')
    if np.any(frames_info['MJD'].values <= date_fix.mjd):
        alt = frames_info['TEL ALT'].values.astype(np.float)
        drot2 = frames_info['INS4 DROT2 BEGIN'].values.astype(np.float)
        pa_correction = np.degrees(np.arctan(np.tan(np.radians(alt-2.*drot2))))
    else:
        pa_correction = 0
    
    # RA/DEC
    ra_drot = frames_info['INS4 DROT2 RA'].values.astype(np.float)
    ra_drot_h = np.floor(ra_drot/1e4)
    ra_drot_m = np.floor((ra_drot - ra_drot_h*1e4)/1e2)
    ra_drot_s = ra_drot - ra_drot_h*1e4 - ra_drot_m*1e2
    ra = coord.Angle((ra_drot_h, ra_drot_m, ra_drot_s), units.hour)
    frames_info['RA'] = ra

    dec_drot = frames_info['INS4 DROT2 DEC'].values.astype(np.float)
    sign = np.sign(dec_drot)
    udec_drot  = np.abs(dec_drot)
    dec_drot_d = np.floor(udec_drot/1e4)
    dec_drot_m = np.floor((udec_drot - dec_drot_d*1e4)/1e2)
    dec_drot_s = udec_drot - dec_drot_d*1e4 - dec_drot_m*1e2
    dec_drot_d *= sign
    dec = coord.Angle((dec_drot_d, dec_drot_m, dec_drot_s), units.degree)
    frames_info['DEC'] = dec
    
    # calculate parallactic angles
    geolon = coord.Angle(frames_info['TEL GEOLON'].values[0], units.degree)
    geolat = coord.Angle(frames_info['TEL GEOLAT'].values[0], units.degree)
    geoelev = frames_info['TEL GEOELEV'].values[0]

    utc = Time(frames_info['TIME START'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG START'] = pa.value + pa_correction

    utc = Time(frames_info['TIME'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG'] = pa.value + pa_correction

    utc = Time(frames_info['TIME END'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG END'] = pa.value + pa_correction

    # pupil offset
    # PA_on-sky = PA_detector + PARANGLE + True_North + PUPOFFSET + IFSOFFSET
    #   PUPOFFSET = 135.99±0.11
    #   IFSOFFSET = 100.48±0.0.10
    drot_mode = frames_info['INS4 DROT2 MODE'].unique()
    if len(drot_mode) != 1:
        raise ValueError('Derotator mode has several values in the sequence')
    if drot_mode == 'ELEV':
        pupoff = 135.99 - 100.48
    elif drot_mode == 'SKY':
        pupoff = -100.48 + frames_info['INS4 DROT2 POSANG']
    elif drot_mode == 'STAT':
        pupoff = -100.48
    else:
        raise ValueError('Unknown derotator mode {0}'.format(drot_mode))

    frames_info['PUPIL OFFSET'] = pupoff

    # final derotation value
    frames_info['DEROT ANGLE'] = frames_info['PARANG'] + pupoff
    

def compute_bad_pixel_map(bpm_files, dtype=np.uint8):
    '''
    Compute a combined bad pixel map provided a list of files

    Parameters
    ----------
    bpm_files : list
        List of names for the bpm files

    dtype : data type
        Data type for the final bpm
    
    Returns
    bpm : array_like
        Combined bad pixel map
    '''

    # check that we have files
    if len(bpm_files) == 0:
        raise ValueError('No bad pixel map files provided')
    
    # get shape
    shape = fits.getdata(bpm_files[0]).shape
    
    # star with empty bpm
    bpm = np.zeros((shape[-2], shape[-1]), dtype=np.uint8)
    
    # fill if files are provided
    for f in bpm_files:
        data = fits.getdata(f)
        bpm = np.logical_or(bpm, data)

    bpm = bpm.astype(dtype)

    return bpm


def collapse_frames_info(finfo, fname, collapse_type, coadd_value=2):
    '''
    Collapse frame info to match the collapse operated on the data

    Parameters
    ----------
    finfo : dataframe
        The data frame with all the information on science frames

    fname : str
       The name of the current file
    
    collapse_type : str
        Type of collapse. Possible values are mean or coadd. Default
        is mean.

    coadd_value : int
        Number of consecutive frames to be coadded when collapse_type
        is coadd. Default is 2

    Returns
    -------
    nfinfo : dataframe
        Collapsed data frame
    '''
    
    print('   ==> collapse frames information')

    nfinfo = None
    if collapse_type == 'none':
        nfinfo = finfo
    elif collapse_type == 'mean':
        index = pd.MultiIndex.from_arrays([[fname], [0]], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index)

        # get min/max indices
        imin = finfo.index.get_level_values(1).min()
        imax = finfo.index.get_level_values(1).max()

        # copy data
        nfinfo.loc[(fname, 0)] = finfo.loc[(fname, imin)]
        
        # update time values
        nfinfo.loc[(fname, 0), 'DET NDIT'] = 1
        nfinfo.loc[(fname, 0), 'TIME START'] = finfo.loc[(fname, imin), 'TIME START']
        nfinfo.loc[(fname, 0), 'TIME END'] = finfo.loc[(fname, imax), 'TIME END']
        nfinfo.loc[(fname, 0), 'TIME'] = finfo.loc[(fname, imin), 'TIME START'] + \
                                         (finfo.loc[(fname, imax), 'TIME END'] - finfo.loc[(fname, imin), 'TIME START']) / 2
        
        # recompute angles
        compute_angles(nfinfo)
    elif collapse_type == 'coadd':
        coadd_value = int(coadd_value)
        NDIT = len(finfo)
        NDIT_new = NDIT // coadd_value

        index = pd.MultiIndex.from_arrays([np.full(NDIT_new, fname), np.arange(NDIT_new)], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index)

        for f in range(NDIT_new):
            # get min/max indices
            imin = int(f*coadd_value)
            imax = int((f+1)*coadd_value-1)

            # copy data
            nfinfo.loc[(fname, f)] = finfo.loc[(fname, imin)]

            # update time values
            nfinfo.loc[(fname, f), 'DET NDIT'] = 1
            nfinfo.loc[(fname, f), 'TIME START'] = finfo.loc[(fname, imin), 'TIME START']
            nfinfo.loc[(fname, f), 'TIME END'] = finfo.loc[(fname, imax), 'TIME END']
            nfinfo.loc[(fname, f), 'TIME'] = finfo.loc[(fname, imin), 'TIME START'] + \
                                             (finfo.loc[(fname, imax), 'TIME END'] - finfo.loc[(fname, imin), 'TIME START']) / 2

        # recompute angles
        compute_angles(nfinfo)
    else:
        raise ValueError('Unknown collapse type {0}'.format(collapse_type))        

    return nfinfo
    

def lines_intersect(a1, a2, b1, b2):
    '''
    Determines the intersection point of two lines passing by points
    (a1,a2) and (b1,b2).

    See https://stackoverflow.com/questions/3252194/numpy-and-line-intersections
    
    Parameters
    ----------
    
    a, b : 2D tuples
        Coordinates of points on line 1
    
    c, d : 2D tuples
        Coordinates of points on line 2
    
    Returns
    -------
    val
        Returns None is lines are parallel, (cx,cy) otherwise.
    '''

    # make sure we have arrays
    a1 = np.array(a1)
    a2 = np.array(a2)
    b1 = np.array(b1)
    b2 = np.array(b2)
    
    # test lines
    da = a2 - a1                # vector from A1 to A2
    db = b2 - b1                # vector from B1 to B2
    dp = a1 - b1
    pda = [-da[1], da[0]]       # perpendicular to A1-A2 vector

    # parallel lines 
    if (pda*db).sum() == 0:
        return None

    # find intersection
    denom = pda @ db
    num   = pda @ dp

    return (num / denom)*db + b1
    

def star_centers_from_waffle_cube(cube, wave, waffle_orientation, high_pass=False, display=False, save_path=None):
    '''
    Compute star center from waffle images

    Parameters
    ----------
    cube : array_like
        Waffle IFS cube

    wave : array_like
        Wavelength values, in nanometers

    waffle_orientation : str
        String giving the waffle orientation '+' or 'x'

    high_pass : bool
        Apply high-pass filter to the image before searching for the satelitte spots
    
    display : bool
        Display the fit of the satelitte spots

    save_path : str
        Path where to save the fit images
    
    Returns
    -------
    spot_center : array_like
        Centers of each individual spot in each frame of the cube

    spot_dist : array_like
        The 6 possible distances between the different spots

    img_center : array_like
        The star center in each frame of the cube
    '''
    # standard parameters
    dim = cube.shape[-1]
    loD = wave*1e-6/8 * 180/np.pi * 3600*1000/pixel
    
    # waffle parameters
    freq = 10 * np.sqrt(2) * 0.97
    box = 8
    if waffle_orientation == '+':
        orient = 0
    elif waffle_orientation == 'x':
        orient = np.pi / 4

    # spot fitting
    xx, yy = np.meshgrid(np.arange(2*box), np.arange(2*box))

    # multi-page PDF to save result
    if save_path is not None:
        pdf = PdfPages(save_path)
    
    # loop over images
    spot_center = np.zeros((nwave, 4, 2))
    spot_dist = np.zeros((nwave, 6))
    ird_center = np.array(((482, 516), (486, 508)))
    img_center = np.full((nwave, 2), ((dim // 2)-1., (dim // 2)-1.))
    for idx, (wave, img) in enumerate(zip(wave, cube)):
        print('  wave {0:2d}/{1:2d} ({2:.3f} micron)'.format(idx+1, nwave, wave))

        # center guess
        cx_int = int(ird_center[idx, 0])
        cy_int = int(ird_center[idx, 1])

        # optional high-pass filter
        if high_pass:
            img = img - ndimage.median_filter(img, 15, mode='mirror')

        # create plot if needed
        if save_path or display:
            fig = plt.figure(0, figsize=(8, 8))
            plt.clf()
            colors = ['red', 'blue', 'magenta', 'purple']
            ax = fig.add_subplot(111)
            ax.imshow(img, aspect='equal', vmin=0, vmax=img.max())
            ax.set_title(r'Image #{0} - {1:.3f} $\mu$m'.format(idx+1, wave))
            
        # satelitte spots
        for s in range(4):
            cx = int(cx_int + freq*loD[idx] * np.cos(orient + np.pi/2*s))
            cy = int(cy_int + freq*loD[idx] * np.sin(orient + np.pi/2*s))

            sub = img[cy-box:cy+box, cx-box:cx+box]

            # fit: Gaussian + constant
            imax = np.unravel_index(np.argmax(sub), sub.shape)
            g_init = models.Gaussian2D(amplitude=sub.max(), x_mean=imax[1], y_mean=imax[0],
                                       x_stddev=loD[idx], y_stddev=loD[idx]) + \
                                       models.Const2D(amplitude=sub.min())
            fitter = fitting.LevMarLSQFitter()
            par = fitter(g_init, xx, yy, sub)
            fit = par(xx, yy)

            cx_final = cx - box + par[0].x_mean
            cy_final = cy - box + par[0].y_mean

            spot_center[idx, s, 0] = cx_final
            spot_center[idx, s, 1] = cy_final

            # plot sattelite spots and fit
            if save_path or display:
                ax.plot([cx_final], [cy_final], marker='D', color=colors[s])
                ax.add_patch(patches.Rectangle((cx-box, cy-box), 2*box, 2*box, ec='white', fc='none'))
                
                axs = fig.add_axes((0.17+s*0.2, 0.17, 0.1, 0.1))
                axs.imshow(sub, aspect='equal', vmin=0, vmax=sub.max())
                axs.plot([par[0].x_mean], [par[0].y_mean], marker='D', color=colors[s])
                axs.set_xticks([])
                axs.set_yticks([])

                axs = fig.add_axes((0.17+s*0.2, 0.06, 0.1, 0.1))
                axs.imshow(fit, aspect='equal', vmin=0, vmax=sub.max())
                axs.set_xticks([])
                axs.set_yticks([])

        # lines intersection
        intersect = lines_intersect(spot_center[idx, 0, :], spot_center[idx, 2, :],
                                    spot_center[idx, 1, :], spot_center[idx, 3, :])
        img_center[idx] = intersect
        
        # scaling
        spot_dist[idx, 0] = np.sqrt(np.sum((spot_center[idx, 0, :] - spot_center[idx, 2, :])**2))
        spot_dist[idx, 1] = np.sqrt(np.sum((spot_center[idx, 1, :] - spot_center[idx, 3, :])**2))
        spot_dist[idx, 2] = np.sqrt(np.sum((spot_center[idx, 0, :] - spot_center[idx, 1, :])**2))
        spot_dist[idx, 3] = np.sqrt(np.sum((spot_center[idx, 0, :] - spot_center[idx, 3, :])**2))
        spot_dist[idx, 4] = np.sqrt(np.sum((spot_center[idx, 1, :] - spot_center[idx, 2, :])**2))
        spot_dist[idx, 5] = np.sqrt(np.sum((spot_center[idx, 2, :] - spot_center[idx, 3, :])**2))

        # finalize plot
        if save_path or display:
            ax.plot([spot_center[idx, 0, 0], spot_center[idx, 2, 0]],
                    [spot_center[idx, 0, 1], spot_center[idx, 2, 1]],
                    color='w', linestyle='dashed')
            ax.plot([spot_center[idx, 1, 0], spot_center[idx, 3, 0]],
                    [spot_center[idx, 1, 1], spot_center[idx, 3, 1]],
                    color='w', linestyle='dashed')

            ax.plot([intersect[0]], [intersect[1]], marker='+', color='w', ms=15)

            ext = 1000 / pixel
            ax.set_xlim(intersect[0]-ext, intersect[0]+ext)
            ax.set_ylim(intersect[1]-ext, intersect[1]+ext)
            
            plt.tight_layout()

            if save_path:
                pdf.savefig()

            if display:
                plt.pause(1e-3)

    if save_path:
        pdf.close()

    return spot_center, spot_dist, img_center


def star_centers_from_PSF_cube(cube, wave, display=False, save_path=None):
    '''
    Compute star center from PSF images

    Parameters
    ----------
    cube : array_like
        PSF IFS cube

    wave : array_like
        Wavelength values, in nanometers

    display : bool
        Display the fit of the satelitte spots

    save_path : str
        Path where to save the fit images
    
    Returns
    -------
    img_center : array_like
        The star center in each frame of the cube
    '''
    
    # standard parameters
    loD = wave*1e-6/8 * 180/np.pi * 3600*1000/pixel
    box = 30
    
    # spot fitting
    xx, yy = np.meshgrid(np.arange(2*box), np.arange(2*box))

    # multi-page PDF to save result
    if save_path is not None:
        pdf = PdfPages(save_path)
        
    # loop over images
    img_center = np.zeros((nwave, 2))
    for idx, (wave, img) in enumerate(zip(wave, cube)):
        print('  wave {0:2d}/{1:2d} ({2:.3f} micron)'.format(idx+1, nwave, wave))

        # center guess
        cy, cx = np.unravel_index(np.argmax(img), img.shape)

        # sub-image
        sub = img[cy-box:cy+box, cx-box:cx+box]

        # fit peak with Gaussian + constant
        imax = np.unravel_index(np.argmax(sub), sub.shape)
        g_init = models.Gaussian2D(amplitude=sub.max(), x_mean=imax[1], y_mean=imax[0],
                                   x_stddev=loD[idx], y_stddev=loD[idx]) + \
                                   models.Const2D(amplitude=sub.min())
        fitter = fitting.LevMarLSQFitter()
        par = fitter(g_init, xx, yy, sub)

        cx_final = cx - box + par[0].x_mean
        cy_final = cy - box + par[0].y_mean

        img_center[idx, 0] = cx_final
        img_center[idx, 1] = cy_final
        
        if save_path or display:
            fig = plt.figure(0, figsize=(8, 8))
            plt.clf()
            ax = fig.add_subplot(111)
            
            ax.imshow(img/img.max(), aspect='equal', vmin=1e-6, vmax=1, norm=colors.LogNorm())
            ax.plot([cx_final], [cy_final], marker='D', color='red')
            ax.add_patch(patches.Rectangle((cx-box, cy-box), 2*box, 2*box, ec='white', fc='none'))
            ax.set_title(r'Image #{0} - {1:.3f} $\mu$m'.format(idx+1, wave))

            ext = 1000 / pixel
            ax.set_xlim(cx_final-ext, cx_final+ext)
            ax.set_ylim(cy_final-ext, cy_final+ext)
                        
            plt.tight_layout()

            if save_path:
                pdf.savefig()

            if display:
                plt.pause(1e-3)

    if save_path:
        pdf.close()

    return img_center


class ImagingReduction(object):
    '''
    SPHERE/IRDIS imaging reduction object
    '''

    ##################################################
    # Class variables
    ##################################################

    
    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, path):
        '''
        Initialization of the IFSReduction

        Parameters
        ----------
        path : str
            Path to the directory containing the raw data
        '''

        # init path
        self._path = ReductionPath.Path(path)

        # execution of recipes
        self._recipe_execution = {
            'sort_files': False,
            'sort_frames': False,
            'check_files_association': False
        }
        
        # reload any existing data frames
        self.read_info()
    
    ##################################################
    # Properties
    ##################################################
    
    @property
    def path(self):
        return self._path

    @property
    def files_info(self):
        return self._files_info
    
    @property
    def frames_info(self):
        return self._frames_info
    
    @property
    def frames_info_preproc(self):
        return self._frames_info_preproc

    @property
    def recipe_execution(self):
        return self._recipe_execution
    
    ##################################################
    # Generic class methods
    ##################################################

    def init_reduction(self):
        '''
        Sort files and frames, perform sanity check
        '''

        # sort files and frames
        self.sort_files()
        self.sort_frames()

        # sanity check
        self.check_files_association()
        
    
    def create_static_calibrations(self):
        '''
        Create all static calibrations with esorex
        '''
        
        self.sph_ird_cal_dark(silent=True)
        self.sph_ird_cal_detector_flat(silent=True)

    
    def preprocess_science(self):
        '''
        Extract images in data cubes
        '''
        
        self.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True,
                                        collapse_science=False, collapse_type='mean', coadd_value=2,
                                        collapse_psf=True, collapse_center=True)


    def process_science(self):
        '''
        Perform star center and combine cubes into final (x,y,time,lambda)
        cubes
        '''
        
        self.sph_ird_star_center(high_pass=False, display=False, save=True)
        self.sph_ird_combine_data(cpix=True, psf_dim=100, science_dim=800, save_scaled=False)

    
    def clean(self):
        '''
        Clean the reduction directory
        '''
        
        self.sph_ird_clean(delete_raw=False, delete_products=False)
        
        
    def full_reduction(self):
        '''
        Performs a full reduction of a data set, from the static
        calibrations to the final (x,y,time,lambda) cubes
        '''
        
        self.init_reduction()
        self.create_static_calibrations()
        self.preprocess_science()
        self.process_science()

    ##################################################
    # SPHERE/IRDIS methods
    ##################################################
    
    def read_info(self):
        '''
        Read the files, calibs and frames information

        files_info : dataframe
            The data frame with all the information on files

        frames_info : dataframe
            The data frame with all the information on science frames

        frames_info_preproc : dataframe
            The data frame with all the information on science frames after pre-processing
        '''

        # path
        path = self._path
        
        # files info
        fname = os.path.join(path.preproc, 'files.csv')
        if os.path.exists(fname):
            files_info = pd.read_csv(fname, index_col=0)

            # convert times
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)

            # update recipe execution
            self._recipe_execution['sort_files'] = True
            if np.any(files_info['PRO CATG'] == 'IRD_MASTER_DARK'):
                self._recipe_execution['sph_ird_cal_dark'] = True
            if np.any(files_info['PRO CATG'] == 'IRD_FLAT_FIELD'):
                self._recipe_execution['sph_ird_cal_detector_flat'] = True
        else:
            files_info = None

        fname = os.path.join(path.preproc, 'frames.csv')
        if os.path.exists(fname):
            frames_info = pd.read_csv(fname, index_col=(0, 1))

            # convert times
            frames_info['DATE-OBS'] = pd.to_datetime(frames_info['DATE-OBS'], utc=True)
            frames_info['DATE'] = pd.to_datetime(frames_info['DATE'], utc=True)
            frames_info['DET FRAM UTC'] = pd.to_datetime(frames_info['DET FRAM UTC'], utc=True)
            frames_info['TIME START'] = pd.to_datetime(frames_info['TIME START'], utc=True)
            frames_info['TIME'] = pd.to_datetime(frames_info['TIME'], utc=True)
            frames_info['TIME END'] = pd.to_datetime(frames_info['TIME END'], utc=True)

            # update recipe execution
            self._recipe_execution['sort_frames'] = True
        else:
            frames_info = None

        fname = os.path.join(path.preproc, 'frames_preproc.csv')
        if os.path.exists(fname):
            frames_info_preproc = pd.read_csv(fname, index_col=(0, 1))

            # convert times
            frames_info_preproc['DATE-OBS'] = pd.to_datetime(frames_info_preproc['DATE-OBS'], utc=True)
            frames_info_preproc['DATE'] = pd.to_datetime(frames_info_preproc['DATE'], utc=True)
            frames_info_preproc['DET FRAM UTC'] = pd.to_datetime(frames_info_preproc['DET FRAM UTC'], utc=True)
            frames_info_preproc['TIME START'] = pd.to_datetime(frames_info_preproc['TIME START'], utc=True)
            frames_info_preproc['TIME'] = pd.to_datetime(frames_info_preproc['TIME'], utc=True)
            frames_info_preproc['TIME END'] = pd.to_datetime(frames_info_preproc['TIME END'], utc=True)            
        else:
            frames_info_preproc = None

        # save data frames in instance variables
        self._files_info = files_info
        self._frames_info = frames_info
        self._frames_info_preproc = frames_info_preproc

        # additional checks to update recipe execution
        if frames_info_preproc is not None:
            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ird_preprocess_science'] = done

            done = True
            files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
                                        (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc_centers'.format(file, idx)
                file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ird_star_center'] = done

        
    def sort_files(self):
        '''
        Sort all raw files and save result in a data frame

        files_info : dataframe
            Data frame with the information on raw files
        '''

        print('Sorting raw files')

        # parameters
        path = self._path
        
        # list files
        files = glob.glob(os.path.join(path.raw, '*.fits'))
        files = [os.path.splitext(os.path.basename(f))[0] for f in files]

        if len(files) == 0:
            raise ValueError('No raw FITS files in root_path directory')

        print(' * found {0} FITS files in {1}'.format(len(files), path.raw))

        # files table
        files_info = pd.DataFrame(index=pd.Index(files, name='FILE'), columns=keywords_short, dtype='float')

        for f in files:
            hdu = fits.open(os.path.join(path.raw, f+'.fits'))
            hdr = hdu[0].header

            for k, sk in zip(keywords, keywords_short):
                files_info.loc[f, sk] = hdr.get(k)

            hdu.close()

        # processed column
        files_info.insert(len(files_info.columns), 'PROCESSED', False)
        files_info.insert(len(files_info.columns), 'PRO CATG', False)

        # convert times
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)

        # save files_info
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))    
        self._files_info = files_info

        # update recipe execution
        self._recipe_execution['sort_files'] = True

        
    def sort_frames(self):
        '''
        Extract the frames information from the science files

        Parameters
        ----------
        root_path : str
            Path to the dataset

        files_info : dataframe
            The data frame with all the information on files

        Returns
        -------
        calibs : dataframe
            A data frame with the information on all frames
        '''

        print('Extracting frames information')

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sort_frames')
        
        # parameters
        path = self._path
        files_info = self._files_info
        
        # science files
        sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] != 'SKY')]    

        # build indices
        files = []
        img   = []
        for file, finfo in sci_files.iterrows():
            NDIT = int(finfo['DET NDIT'])

            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))

        # create new dataframe
        frames_info = pd.DataFrame(columns=sci_files.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        frames_info = frames_info.align(files_info, level=0)[1]    

        # compute timestamps
        compute_times(frames_info)

        # compute angles (ra, dec, parang)
        compute_angles(frames_info)

        # save
        frames_info.to_csv(os.path.join(path.preproc, 'frames.csv'))
        self._frames_info = frames_info

        # update recipe execution
        self._recipe_execution['sort_frames'] = True


    def check_files_association(self):
        '''
        Performs the calibration files association as a sanity check
        '''

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'check_files_association')
        
        print('Performing file association for calibrations')

        # parameters
        files_info = self._files_info

        # instrument arm
        arm = files_info['SEQ ARM'].unique()
        if len(arm) != 1:
            raise ValueError('Sequence is mixing different instruments: {0}'.format(arm))
        
        # IRDIS obs mode and filter combination
        modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'].unique()
        if len(modes != 1):
            raise ValueError('Sequence is mixing different types of observations: {0}'.format(modes))

        filter_combs = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()
        if len(filter_combs != 1):
            raise ValueError('Sequence is mixing different types of filters combinations: {0}'.format(filter_combs))
        filter_comb = filter_combs[0]
        
        # specific data frame for calibrations
        # keep static calibrations and sky backgrounds
        calibs = files_info[(files_info['DPR CATG'] == 'CALIB') |
                            ((files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == 'SKY'))]

        ###############################################
        # static calibrations not dependent on science
        ###############################################
        error_flag = 0
        warning_flag = 0

        # flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS COMB IFLT'] == filter_comb)]
        if len(cfiles) <= 1:
            error_flag += 1
            print(' * Error: there should be more than 1 flat in filter combination {0}'.format(filter_comb))

        ##################################################
        # static calibrations that depend on science (DIT)
        ##################################################

        obj = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'DPR TYPE'].apply(lambda s: s[0:6])
        DITs = files_info.loc[(files_info['DPR CATG'] == 'SCIENCE') & (obj == 'OBJECT'), 'DET SEQ1 DIT'].unique().round(2)

        # handle darks in a slightly different way because there might be several different DITs
        for DIT in DITs:
            # instrumental backgrounds
            cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                            (calibs['DET SEQ1 DIT'].round(2) == DIT)]
            if len(cfiles) == 0:
                warning_flag += 1
                print(' * Warning: there is no dark/background for science files with DIT={0} sec. '.format(DIT) +
                      'it is *highly recommended* to include one to obtain the best data reduction. ' +
                      'A single dark/background file is sufficient, and it can easily be downloaded ' +
                      'from the ESO archive')

            # sky backgrounds
            cfiles = files_info[(files_info['DPR TYPE'] == 'SKY') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
            if len(cfiles) == 0:
                warning_flag += 1
                print(' * Warning: there is no sky background for science files with DIT={0} sec. '.format(DIT) +
                      'Using a sky background instead of an internal instrumental background can ' +
                      'usually provide a cleaner data reduction, especially in K-band')

        # error reporting
        print('There are {0} warning(s) and {1} error(s) in the classification of files'.format(warning_flag, error_flag))
        if error_flag:
            raise ValueError('There is {0} errors that should be solved before proceeding'.format(error_flag))

        
    def sph_ird_cal_dark(self, silent=True):
        '''
        Create the dark and background calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Optional, default is True
        '''

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sph_ird_cal_dark')
        
        print('Creating darks and backgrounds')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        raw = files_info[np.logical_not(files_info['PROCESSED'])]
        calibs = raw[(files_info['DPR TYPE'] == 'DARK') | (files_info['DPR TYPE'] == 'DARK,BACKGROUND') |
                     (files_info['DPR TYPE'] == 'SKY')]

        # loops on type and DIT value
        types = ['DARK', 'DARK,BACKGROUND', 'SKY']
        DITs = calibs['DET SEQ1 DIT'].unique().round(2)
        filter_combs = calibs['INS COMB IFLT'].unique()

        for ctype in types:
            for DIT in DITs:
                for cfilt in filter_combs:
                    cfiles = calibs[(calibs['DPR TYPE'] == ctype) & (calibs['DET SEQ1 DIT'].round(2) == DIT) &
                                    (calibs['INS COMB IFLT'] == cfilt)]
                    files = cfiles.index

                    # skip non-existing combinations
                    if len(cfiles) == 0:
                        continue

                    print(' * {0} in filter {1} with DIT={2:.2f} sec ({3} files)'.format(ctype, cfilt, DIT, len(cfiles)))

                    # create sof
                    sof = os.path.join(path.sof, 'dark_filt={0}_DIT={1:.2f}.sof'.format(cfilt, DIT))
                    file = open(sof, 'w')
                    for f in files:
                        file.write('{0}{1}.fits     {2}\n'.format(path.raw, f, 'IRD_DARK_RAW'))
                    file.close()

                    # products
                    if ctype == 'SKY':
                        loc = 'sky'
                    else:
                        loc = 'internal'
                    dark_file = 'dark_{0}_filt={1}_DIT={2:.2f}'.format(loc, cfilt, DIT)
                    bpm_file  = 'dark_{0}_bpm_filt={1}_DIT={2:.2f}'.format(loc, cfilt, DIT)

                    # different max level in K-band
                    max_level = 1000
                    if cfilt in ['DB_K12', 'BB_Ks']:
                        max_level = 15000
                    
                    # esorex parameters    
                    args = ['esorex',
                            '--no-checksum=TRUE',
                            '--no-datamd5=TRUE',
                            'sph_ird_master_dark',
                            '--ird.master_dark.sigma_clip=5.0',
                            '--ird.master_dark.save_addprod=TRUE',
                            '--ird.master_dark.max_acceptable={0}'.format(max_level),
                            '--ird.master_dark.outfilename={0}{1}.fits'.format(path.calib, dark_file),
                            '--ird.master_dark.badpixfilename={0}{1}.fits'.format(path.calib, bpm_file),
                            sof]

                    # check esorex
                    if shutil.which('esorex') is None:
                        raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                                        'that the ESO pipeline is properly installed before running pySPHERE.')

                    # execute esorex
                    if silent:
                        proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
                    else:
                        proc = subprocess.run(args, cwd=path.tmp)

                    if proc.returncode != 0:
                        raise ValueError('esorex process was not successful')

                    # store products
                    files_info.loc[dark_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
                    files_info.loc[dark_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
                    files_info.loc[dark_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'][0]
                    files_info.loc[dark_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'][0]
                    files_info.loc[dark_file, 'INS1 MODE'] = cfiles['INS1 MODE'][0]
                    files_info.loc[dark_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'][0]
                    files_info.loc[dark_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'][0]
                    files_info.loc[dark_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
                    files_info.loc[dark_file, 'PROCESSED'] = True
                    files_info.loc[dark_file, 'PRO CATG'] = 'IRD_MASTER_DARK'

                    files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
                    files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
                    files_info.loc[bpm_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'][0]
                    files_info.loc[bpm_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'][0]
                    files_info.loc[bpm_file, 'INS1 MODE'] = cfiles['INS1 MODE'][0]
                    files_info.loc[bpm_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'][0]
                    files_info.loc[bpm_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'][0]
                    files_info.loc[bpm_file, 'PROCESSED'] = True
                    files_info.loc[bpm_file, 'PRO CATG']  = 'IRD_STATIC_BADPIXELMAP'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ird_cal_dark'] = True


    def sph_ird_cal_detector_flat(self, silent=True):
        '''
        Create the detector flat calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Optional, default is True
        '''

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sph_ird_cal_detector_flat')
        
        print('Creating flats')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        raw = files_info[np.logical_not(files_info['PROCESSED'])]
        calibs = raw[(files_info['DPR TYPE'] == 'FLAT,LAMP') | (files_info['DPR TECH'] == 'IMAGE')]
        filter_combs = calibs['INS COMB IFLT'].unique()
        
        for cfilt in filter_combs:
            cfiles = calibs[calibs['INS COMB IFLT'] == cfilt]
            files = cfiles.index

            print(' * filter {0} ({1} files)'.format(cfilt, len(cfiles)))
            
            # create sof
            sof = os.path.join(path.sof, 'flat_filt={0}.sof'.format(cfilt))
            file = open(sof, 'w')
            for f in files:
                file.write('{0}{1}.fits     {2}\n'.format(path.raw, f, 'IRD_FLAT_FIELD_RAW'))
            file.close()

            # products
            flat_file = 'flat_filt={0}'.format(cfilt)
            bpm_file  = 'flat_bpm_filt={0}'.format(cfilt)
            
            # esorex parameters    
            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_instrument_flat',
                    '--ird.instrument_flat.save_addprod=TRUE',
                    '--ird.instrument_flat.outfilename={0}{1}.fits'.format(path.calib, flat_file),
                    '--ird.instrument_flat.badpixfilename={0}{1}.fits'.format(path.calib, bpm_file),
                    sof]

            # check esorex
            if shutil.which('esorex') is None:
                raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                                'that the ESO pipeline is properly installed before running pySPHERE.')

            # execute esorex
            if silent:
                proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
            else:
                proc = subprocess.run(args, cwd=path.tmp)

            if proc.returncode != 0:
                raise ValueError('esorex process was not successful')

            # store products
            files_info.loc[flat_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            files_info.loc[flat_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            files_info.loc[flat_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'][0]
            files_info.loc[flat_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'][0]
            files_info.loc[flat_file, 'INS1 MODE'] = cfiles['INS1 MODE'][0]
            files_info.loc[flat_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'][0]
            files_info.loc[flat_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'][0]
            files_info.loc[flat_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
            files_info.loc[flat_file, 'PROCESSED'] = True
            files_info.loc[flat_file, 'PRO CATG'] = 'IRD_FLAT_FIELD'

            files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            files_info.loc[bpm_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'][0]
            files_info.loc[bpm_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'][0]
            files_info.loc[bpm_file, 'INS1 MODE'] = cfiles['INS1 MODE'][0]
            files_info.loc[bpm_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'][0]
            files_info.loc[bpm_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'][0]
            files_info.loc[bpm_file, 'PROCESSED'] = True
            files_info.loc[bpm_file, 'PRO CATG']  = 'IRD_NON_LINEAR_BADPIXELMAP'
        
        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_detector_flat'] = True


    def sph_ird_preprocess_science(self,
                                   subtract_background=True, fix_badpix=True,
                                   collapse_science=False, collapse_type='mean', coadd_value=2,
                                   collapse_psf=True, collapse_center=True):
        '''
        Pre-processes the science frames.

        This function can perform multiple steps:
          - collapse of the frames according to different schemes
          - subtract the background
          - correct bad pixels
          - reformat IRDIS data in (x,y,lambda) cubes

        For the science, 2 collapse methods are available: mean or
        coadd. With mean, the full cubes are mean-combined into a single
        frame. With coadd, the frames are coadded following the
        coadd_value. This can result in lost frames if the number of NDIT
        is not a multiple of coadd_value.

        For the PSFs and star center frames, there is either no collapse
        or a mean collapse.

        Parameters
        ----------
        subtract_background : bool
            Performs background subtraction. Default is True

        fix_badpix : bool
            Performs correction of bad pixels. Default is True

        collapse_science :  bool
            Collapse data for OBJECT cubes. Default is False

        collapse_type : str
            Type of collapse. Possible values are mean or coadd. Default
            is mean.

        coadd_value : int
            Number of consecutive frames to be coadded when collapse_type
            is coadd. Default is 2

        collapse_psf :  bool
            Collapse data for OBJECT,FLUX cubes. Default is True. Note
            that the collapse type is mean and cannot be changed.

        collapse_center :  bool
            Collapse data for OBJECT,CENTER cubes. Default is True. Note
            that the collapse type is mean and cannot be changed.    
        '''

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sph_ird_preprocess_science')
        
        print('Pre-processing science files')

        # parameters
        path = self._path
        files_info = self._files_info
        frames_info = self._frames_info
        
        # clean before we start
        files = glob.glob(os.path.join(path.preproc, '*_DIT???_preproc.fits'))
        for file in files:
            os.remove(file)

        # filter combination
        filter_comb = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()[0]

        # bpm
        if fix_badpix:
            bpm_files = files_info[(files_info['PRO CATG'] == 'IRD_STATIC_BADPIXELMAP') |
                                   (files_info['PRO CATG'] == 'IRD_NON_LINEAR_BADPIXELMAP')].index
            bpm_files = [os.path.join(path.calib, f+'.fits') for f in bpm_files]

            bpm = compute_bad_pixel_map(bpm_files)

        # flat        
        flat_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_FLAT_FIELD') &
                               (files_info['INS COMB IFLT'] == filter_comb)]
        if len(flat_file) != 1:
            raise ValueError('There should be exactly 1 flat file. Found {0}.'.format(len(flat_file)))
        flat = fits.getdata(os.path.join(path.calib, flat_file.index[0]+'.fits'))
            
        # final dataframe
        index = pd.MultiIndex(names=['FILE', 'IMG'], levels=[[], []], labels=[[], []])
        frames_info_preproc = pd.DataFrame(index=index, columns=frames_info.columns)

        # loop on the different type of science files
        sci_types = ['OBJECT,CENTER', 'OBJECT,FLUX', 'OBJECT']
        dark_types = ['SKY', 'DARK,BACKGROUND', 'DARK']
        for typ in sci_types:
            # science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == typ)]
            sci_DITs = list(sci_files['DET SEQ1 DIT'].round(2).unique())

            if len(sci_files) == 0:
                continue        

            for DIT in sci_DITs:
                sfiles = sci_files[sci_files['DET SEQ1 DIT'].round(2) == DIT]

                print('{0} files of type {1} with DIT={2} sec'.format(len(sfiles), typ, DIT))

                if subtract_background:
                    # look for sky, then background, then darks
                    # normally there should be only one with the proper DIT
                    dfiles = []
                    for d in dark_types:
                        dfiles = files_info[(files_info['PRO CATG'] == 'IRD_MASTER_DARK') &
                                            (files_info['DPR TYPE'] == d) & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
                        if len(dfiles) != 0:
                            break
                    print('   ==> found {0} corresponding {1} file'.format(len(dfiles), d))

                    if len(dfiles) != 1:
                        raise ValueError('Unexpected number of background files ({0})'.format(len(dfiles)))

                    bkg = fits.getdata(os.path.join(path.calib, dfiles.index[0]+'.fits'))

                # process files
                for idx, (fname, finfo) in enumerate(sci_files.iterrows()):
                    # frames_info extract
                    finfo = frames_info.loc[(fname, slice(None)), :]

                    print(' * file {0}/{1}: {2}, NDIT={3}'.format(idx+1, len(sci_files), fname, len(finfo)))

                    # read data
                    print('   ==> read data')
                    img, hdr = fits.getdata(os.path.join(path.raw, fname+'.fits'), header=True)
                    
                    # add extra dimension to single images to make cubes
                    if img.ndim == 2:
                        img = img[np.newaxis, ...]
                        
                    # collapse
                    if (typ == 'OBJECT,CENTER'):
                        if collapse_center:
                            print('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                        else:
                            frames_info_new = collapse_frames_info(finfo, fname, 'none')
                    elif (typ == 'OBJECT,FLUX'):
                        if collapse_psf:
                            print('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                        else:
                            frames_info_new = collapse_frames_info(finfo, fname, 'none')
                    elif (typ == 'OBJECT'):
                        if collapse_science:                        
                            if collapse_type == 'mean':
                                print('   ==> collapse: mean ({0} -> 1 frame, 0 dropped)'.format(len(img)))
                                img = np.mean(img, axis=0, keepdims=True)

                                frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                            elif collapse_type == 'coadd':
                                if (not isinstance(coadd_value, int)) or (coadd_value <= 1):
                                    raise TypeError('coadd_value must be an integer >1')

                                coadd_value = int(coadd_value)
                                NDIT = len(img)
                                NDIT_new = NDIT // coadd_value
                                dropped = NDIT % coadd_value

                                if coadd_value > NDIT:
                                    raise ValueError('coadd_value ({0}) must be < NDIT ({1})'.format(coadd_value, NDIT))

                                print('   ==> collapse: mean ({0} -> {1} frames, {2} dropped)'.format(NDIT, NDIT_new, dropped))

                                # coadd frames
                                nimg = np.empty((NDIT_new, 2048, 2048), dtype=img.dtype)
                                for f in range(NDIT_new):
                                    nimg[f] = np.mean(img[f*coadd_value:(f+1)*coadd_value], axis=0)
                                img = nimg

                                frames_info_new = collapse_frames_info(finfo, fname, 'coadd', coadd_value=coadd_value)
                            else:
                                raise ValueError('Unknown collapse type {0}'.format(collapse_type))
                        else:
                            frames_info_new = collapse_frames_info(finfo, fname, 'none')

                    # merge collapse collapsed frames_info
                    frames_info_preproc = pd.concat((frames_info_preproc, frames_info_new))

                    # background subtraction
                    if subtract_background:
                        print('   ==> subtract background')
                        for f in range(len(img)):
                            img[f] -= bkg

                    # divide flat
                    if subtract_background:
                        print('   ==> divide by flat field')
                        for f in range(len(img)):
                            img[f] /= flat
                        
                    # bad pixels correction
                    if fix_badpix:
                        print('   ==> correct bad pixels')
                        for f in range(len(img)):                            
                            frame = img[f]
                            frame = imutils.fix_badpix(frame, bpm)
                            img[f] = frame

                    # reshape data
                    print('   ==> reshape data')
                    NDIT = img.shape[0]
                    nimg = np.zeros((NDIT, 2, 1024, 1024))
                    for f in range(len(img)):
                        nimg[f, 0] = img[f, :, 0:1024]
                        nimg[f, 1] = img[f, :, 1024:]
                    img = nimg
                        
                    # save DITs individually
                    for f in range(len(img)):
                        frame = nimg[f, ...].squeeze()                    
                        hdr['HIERARCH ESO DET NDIT'] = 1
                        fits.writeto(os.path.join(path.preproc, fname+'_DIT{0:03d}_preproc.fits'.format(f)), frame, hdr,
                                     overwrite=True, output_verify='silentfix')

                    print()

            print()

        # sort and save final dataframe
        frames_info_preproc.sort_values(by='TIME', inplace=True)
        frames_info_preproc.to_csv(os.path.join(path.preproc, 'frames_preproc.csv'))

        self._frames_info_preproc = frames_info_preproc

        # update recipe execution
        self._recipe_execution['sph_ird_preprocess_science'] = True


    def sph_ird_star_center(self, high_pass=False, display=False, save=True):
        '''
        Determines the star center for all frames

        Parameters
        ----------
        high_pass : bool
            Apply high-pass filter to the image before searching for the satelitte spots

        display : bool
            Display the fit of the satelitte spots

        save : bool
            Save the fit of the sattelite spot for quality check. Default is True,
            although it is a bit slow.
        '''

        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sph_ird_star_center')
        
        print('Star centers determination')

        # parameters
        path = self._path
        frames_info = self._frames_info_preproc

        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)                
        wave = np.array(wave) / 1000.
        
        # start with OBJECT,FLUX
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']        
        if len(flux_files) != 0:
            for file, idx in flux_files.index:
                print('  ==> OBJECT,FLUX: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube, hdr = fits.getdata(files[0], header=True)

                # centers
                if save:
                    save_path = os.path.join(path.products, fname+'_PSF_fitting.pdf')
                else:
                    save_path = None
                img_center = star_centers_from_PSF_cube(cube, wave, display=display, save_path=save_path)

                # save
                fits.writeto(os.path.join(path.preproc, fname+'_centers.fits'), img_center, overwrite=True)
                print()

        # then OBJECT,CENTER
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) != 0:
            for file, idx in starcen_files.index:
                print('  ==> OBJECT,CENTER: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube, hdr = fits.getdata(files[0], header=True)
                
                # centers
                waffle_orientation = hdr['HIERARCH ESO OCS WAFFLE ORIENT']
                if save:
                    save_path = os.path.join(path.products, fname+'_spots_fitting.pdf')
                else:
                    save_path = None
                spot_center, spot_dist, img_center = star_centers_from_waffle_cube(cube, wave, waffle_orientation,
                                                                                   high_pass=high_pass, display=display,
                                                                                   save_path=save_path)

                # save
                fits.writeto(os.path.join(path.preproc, fname+'_centers.fits'), img_center, overwrite=True)
                print()

        # update recipe execution
        self._recipe_execution['sph_ifs_star_center'] = True


    def sph_ird_combine_data(self, cpix=True, psf_dim=80, science_dim=290, save_scaled=False):
        '''
        Combine and save the science data into final cubes
        
        Parameters
        ----------
        cpix : bool
            If True the images are centered on the pixel at coordinate
            (dim//2,dim//2). If False the images are centered between 4
            pixels, at coordinates ((dim-1)/2,(dim-1)/2). Default is True.

        psf_dim : even int
            Size of the PSF images. Default is 80x80 pixels

        science_dim : even int    
            Size of the science images (star centers and standard
            coronagraphic images). Default is 290, 290 pixels

        save_scaled : bool    
            Also save the wavelength-rescaled cubes. Makes the process
            much longer. The default is False

        '''
        
        # check if recipe can be executed
        check_recipe_execution(self._recipe_execution, 'sph_ird_combine_data')
        
        print('Combine science data')

        # parameters
        path = self._path
        frames_info = self._frames_info_preproc
        
        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)                
        wave = np.array(wave) / 1000.        

        fits.writeto(os.path.join(path.products, 'wavelength.fits'), wave, overwrite=True)

        # max images size
        if psf_dim > 1024:
            print('Warning: psf_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            psf_dim = 1024

        if science_dim > 1024:
            print('Warning: science_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            science_dim = 1024
            
        #
        # frames info
        #
        frames_info.to_csv(os.path.join(path.products, 'frames.csv'))

        #
        # OBJECT,FLUX
        #
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        nfiles = len(flux_files)
        if nfiles != 0:
            print(' * OBJECT,FLUX data')

            # final arrays
            psf_cube   = np.zeros((nwave, nfiles, psf_dim, psf_dim))
            psf_parang = np.zeros(nfiles)
            psf_derot  = np.zeros(nfiles)
            if save_scaled:
                psf_cube_scaled = np.zeros((nwave, nfiles, psf_dim, psf_dim))

            # final center
            if cpix:
                cc = psf_dim // 2
            else:
                cc = (psf_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(flux_files.index):
                print('  ==> {0}, DIT #{1}'.format(file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube = fits.getdata(files[0])
                centers = fits.getdata(os.path.join(path.preproc, fname+'_centers.fits'))

                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave*1000)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                psf_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                psf_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method='fft')
                    nimg = nimg / DIT / attenuation[wave_idx]

                    psf_cube[wave_idx, file_idx] = nimg[:psf_dim, :psf_dim]

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = psf_cube[wave_idx, file_idx]
                        psf_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

            # save final cubes
            fits.writeto(os.path.join(path.products, 'psf_cube.fits'), psf_cube, overwrite=True)
            fits.writeto(os.path.join(path.products, 'psf_parang.fits'), psf_parang, overwrite=True)
            fits.writeto(os.path.join(path.products, 'psf_derot.fits'), psf_derot, overwrite=True)
            if save_scaled:
                fits.writeto(os.path.join(path.products, 'psf_cube_scaled.fits'), psf_cube_scaled, overwrite=True)

            # delete big cubes
            del psf_cube
            if save_scaled:
                del psf_cube_scaled

            print()

        #
        # OBJECT,CENTER
        #
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        nfiles = len(starcen_files)
        if nfiles != 0:
            print(' * OBJECT,CENTER data')

            # final arrays
            cen_cube   = np.zeros((nwave, nfiles, science_dim, science_dim))
            cen_parang = np.zeros(nfiles)
            cen_derot  = np.zeros(nfiles)
            if save_scaled:
                cen_cube_scaled = np.zeros((nwave, nfiles, science_dim, science_dim))

            # final center
            if cpix:
                cc = science_dim // 2
            else:
                cc = (science_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(starcen_files.index):
                print('  ==> {0}, DIT #{1}'.format(file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube = fits.getdata(files[0])
                centers = fits.getdata(os.path.join(path.preproc, fname+'_centers.fits'))

                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave*1000)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                cen_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                cen_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method='fft')
                    nimg = nimg / DIT / attenuation[wave_idx]

                    cen_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = cen_cube[wave_idx, file_idx]
                        cen_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

            # save final cubes
            fits.writeto(os.path.join(path.products, 'star_center_cube.fits'), cen_cube, overwrite=True)
            fits.writeto(os.path.join(path.products, 'star_center_parang.fits'), cen_parang, overwrite=True)
            fits.writeto(os.path.join(path.products, 'star_center_derot.fits'), cen_derot, overwrite=True)
            if save_scaled:
                fits.writeto(os.path.join(path.products, 'star_center_cube_scaled.fits'), cen_cube_scaled, overwrite=True)

            # delete big cubes
            del cen_cube
            if save_scaled:
                del cen_cube_scaled

            print()

        #
        # OBJECT
        #
        object_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT']
        nfiles = len(object_files)
        if nfiles != 0:
            print(' * OBJECT data')

            # get first DIT of first OBJECT,CENTER in the sequence. See issue #12.
            starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
            if len(starcen_files) == 0:
                print(' ==> no OBJECT,CENTER file in the data set. Images cannot be accurately centred. ' +
                      'They will just be combined.')

                centers = np.full((nwave, 2), cc)
            else:
                fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                centers = fits.getdata(os.path.join(path.preproc, fname))

            # final arrays
            sci_cube   = np.zeros((nwave, nfiles, science_dim, science_dim))
            sci_parang = np.zeros(nfiles)
            sci_derot  = np.zeros(nfiles)
            if save_scaled:
                sci_cube_scaled = np.zeros((nwave, nfiles, science_dim, science_dim))

            # final center
            if cpix:
                cc = science_dim // 2
            else:
                cc = (science_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(object_files.index):
                print('  ==> {0}, DIT #{1}'.format(file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube = fits.getdata(files[0])

                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave*1000)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                sci_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                sci_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method='fft')
                    nimg = nimg / DIT / attenuation[wave_idx]

                    sci_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = sci_cube[wave_idx, file_idx]
                        sci_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

            # save final cubes
            fits.writeto(os.path.join(path.products, 'science_cube.fits'), sci_cube, overwrite=True)
            fits.writeto(os.path.join(path.products, 'science_parang.fits'), sci_parang, overwrite=True)
            fits.writeto(os.path.join(path.products, 'science_derot.fits'), sci_derot, overwrite=True)
            if save_scaled:
                fits.writeto(os.path.join(path.products, 'science_cube_scaled.fits'), sci_cube_scaled, overwrite=True)

            # delete big cubes
            del sci_cube
            if save_scaled:
                del sci_cube_scaled

            print()

        # update recipe execution
        self._recipe_execution['sph_ifs_combine_data'] = True

            

    def sph_ird_clean(self, delete_raw=False, delete_products=False):
        '''
        Clean everything except for raw data and science products (by default)

        Parameters
        ----------
        delete_raw : bool
            Delete raw data. Default is False

        delete_products : bool
            Delete science products. Default is False
        '''

        # parameters
        path = self._path
                
        # tmp
        if os.path.exists(path.tmp):
            shutil.rmtree(path.tmp, ignore_errors=True)

        # sof
        if os.path.exists(path.sof):
            shutil.rmtree(path.sof, ignore_errors=True)

        # calib
        if os.path.exists(path.calib):
            shutil.rmtree(path.calib, ignore_errors=True)

        # preproc
        if os.path.exists(path.preproc):
            shutil.rmtree(path.preproc, ignore_errors=True)

        # raw
        if delete_raw:
            if os.path.exists(path.raw):
                shutil.rmtree(path.raw, ignore_errors=True)

        # products
        if delete_products:
            if os.path.exists(path.products):
                shutil.rmtree(path.products, ignore_errors=True)
