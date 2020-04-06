import pandas as pd
import numpy as np
import astropy.coordinates as coord
import astropy.units as units
import scipy.ndimage as ndimage
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import logging

import sphere
import sphere.utils.aperture as aperture

from astropy.io import fits
from astropy.time import Time
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages

global_cmap = 'inferno'

_log = logging.getLogger(__name__)


def recipe_executable(recipes_status, reduction_status, recipe, requirements, logger=_log):
    '''
    Check if a recipe is executabled given the status of other recipes

    Parameters
    ----------
    recipes_status : dict
        Status of executed recipes

    reduction_status : sphere state
        Overall status of the reduction

    recipe : str
        Name of the current recipe

    requirements : dict
        Dictionary providing the recipe requirements

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    execute_recipe : bool
        Current recipe can be executed safely
    '''
    
    if reduction_status == sphere.FATAL:
        logger.critical('   ==> reduction is in a FATAL state! See log file for details')
        return False
    
    recipes = recipes_status.keys()
    requirements = requirements[recipe]
    
    execute_recipe = True
    missing = []
    for r in requirements:
        if r not in recipes:
            execute_recipe = False
            missing.append(r)
        elif recipes_status[r] != sphere.SUCCESS:
            execute_recipe = False
            missing.append(r)

    if not execute_recipe:
        logger.error('{} cannot be executed because the following recipes have not been executed or have result in unrecoverable errors: {}. '.format(recipe, missing))
        recipes_status[recipe] = sphere.ERROR

    logger.debug('> execution requirements check for {}: {}'.format(recipe, execute_recipe))
    
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


def compute_times(frames_info, logger=_log):
    '''
    Compute the various timestamps associated to frames

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''

    logger.debug('> compute time stamps')
    
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


def compute_angles(frames_info, logger=_log):
    '''
    Compute the various angles associated to frames: RA, DEC, parang,
    pupil offset, final derotation angle

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''

    logger.debug('> compute angles')
    
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
    ra_hour = coord.Angle((ra_drot_h, ra_drot_m, ra_drot_s), units.hour)
    ra_deg  = ra_hour*15
    frames_info['RA'] = ra_deg

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
    ha  = lst - ra_hour
    pa  = parallatic_angle(ha, dec[0], geolat)
    frames_info['PARANG START'] = pa.value + pa_correction
    frames_info['HOUR ANGLE START'] = ha
    frames_info['LST START'] = lst

    utc = Time(frames_info['TIME'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra_hour
    pa  = parallatic_angle(ha, dec[0], geolat)
    frames_info['PARANG'] = pa.value + pa_correction
    frames_info['HOUR ANGLE'] = ha
    frames_info['LST'] = lst

    utc = Time(frames_info['TIME END'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra_hour
    pa  = parallatic_angle(ha, dec[0], geolat)
    frames_info['PARANG END'] = pa.value + pa_correction
    frames_info['HOUR ANGLE END'] = ha
    frames_info['LST END'] = lst

    #
    # Derotation angles
    #
    # PA_on-sky = PA_detector + PARANGLE + True_North + PUP_OFFSET + INSTRUMENT_OFFSET
    #  PUP_OFFSET = -135.99 ± 0.11
    #  INSTRUMENT_OFFSET
    #   IFS = +100.48 ± 0.10
    #   IRD =    0.00 ± 0.00
    #
    instru = frames_info['SEQ ARM'].unique()
    if len(instru) != 1:
        logger.error('Sequence is mixing different instruments: {0}'.format(instru))
        return sphere.ERROR
    if instru == 'IFS':
        instru_offset = -100.48
    elif instru == 'IRDIS':
        instru_offset = 0.0
    else:
        logger.error('Unkown instrument {0}'.format(instru))
        return sphere.ERROR

    drot_mode = frames_info['INS4 DROT2 MODE'].unique()
    if len(drot_mode) != 1:
        logger.error('Derotator mode has several values in the sequence')
        return sphere.ERROR
    if drot_mode == 'ELEV':
        pupoff = 135.99
    elif drot_mode == 'SKY':
        pupoff = -100.48 + frames_info['INS4 DROT2 POSANG']
    elif drot_mode == 'STAT':
        pupoff = -100.48
    else:
        logger.error('Unknown derotator mode {0}'.format(drot_mode))
        return sphere.ERROR

    frames_info['PUPIL OFFSET'] = pupoff + instru_offset

    # final derotation value
    frames_info['DEROT ANGLE'] = frames_info['PARANG'] + pupoff + instru_offset
    
    return sphere.SUCCESS


def compute_bad_pixel_map(bpm_files, dtype=np.uint8, logger=_log):
    '''
    Compute a combined bad pixel map provided a list of files

    Parameters
    ----------
    bpm_files : list
        List of names for the bpm files

    dtype : data type
        Data type for the final bpm

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    bpm : array_like
        Combined bad pixel map
    '''

    logger.debug('> compute master bad pixel map from {} files'.format(len(bpm_files)))
    
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


def collapse_frames_info(finfo, fname, collapse_type, coadd_value=2, logger=_log):
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

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    nfinfo : dataframe
        Collapsed data frame, or None in case of error
    '''

    logger.info('   ==> collapse frames information')

    nfinfo = None
    if collapse_type == 'none':
        nfinfo = finfo
        logger.debug('> type=none: copy input data frame')
    elif collapse_type == 'mean':
        index = pd.MultiIndex.from_arrays([[fname], [0]], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index, dtype=np.float)

        logger.debug('> type=mean: extract min/max values')
        
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
        ret = compute_angles(nfinfo, logger=logger)
        if ret == sphere.ERROR:
            return None
    elif collapse_type == 'coadd':
        coadd_value = int(coadd_value)
        NDIT = len(finfo)
        NDIT_new = NDIT // coadd_value

        logger.debug('> type=coadd: extract sub-groups of {} frames'.format(coadd_value))
        
        index = pd.MultiIndex.from_arrays([np.full(NDIT_new, fname), np.arange(NDIT_new)], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index, dtype=np.float)

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
        ret = compute_angles(nfinfo, logger=logger)
        if ret == sphere.ERROR:
            return None
    else:
        logger.error('Unknown collapse type {0}'.format(collapse_type))
        return None

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


def star_centers_from_PSF_img_cube(cube, wave, pixel, exclude_fraction=0.1,
                                   save_path=None, logger=_log):
    '''
    Compute star center from PSF images (IRDIS CI, IRDIS DBI, IFS)

    Parameters
    ----------
    cube : array_like
        IRDIFS PSF cube

    wave : array_like
        Wavelength values, in nanometers

    pixel : float
        Pixel scale, in mas/pixel

    exclude_fraction : float
        Exclude a fraction of the image borders to avoid getting
        biased by hot pixels close to the edges. Default is 10%

    save_path : str
        Path where to save the fit images. Default is None, which means
        that the plot is not produced

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    img_centers : array_like
        The star center in each frame of the cube

    '''

    # standard parameters
    nwave = wave.size
    loD = wave*1e-9/8 * 180/np.pi * 3600*1000/pixel
    box = 30

    # spot fitting
    xx, yy = np.meshgrid(np.arange(2*box), np.arange(2*box))

    # loop over images
    img_centers    = np.zeros((nwave, 2))
    failed_centers = np.zeros(nwave, dtype=np.bool)
    for idx, (cwave, img) in enumerate(zip(wave, cube)):
        logger.info('   ==> wave {0:2d}/{1:2d} ({2:4.0f} nm)'.format(idx+1, nwave, cwave))

        # remove any NaN
        img = np.nan_to_num(img)        
        
        # center guess
        cy, cx = np.unravel_index(np.argmax(img), img.shape)

        # check if we are really too close to the edge
        dim = img.shape
        lf = exclude_fraction
        hf = 1-exclude_fraction
        if (cx <= lf*dim[-1]) or (cx >= hf*dim[-1]) or \
           (cy <= lf*dim[0])  or (cy >= hf*dim[0]):
            nimg = img.copy()
            nimg[:, :int(lf*dim[-1])] = 0
            nimg[:, int(hf*dim[-1]):] = 0
            nimg[:int(lf*dim[0]), :]  = 0
            nimg[int(hf*dim[0]):, :]  = 0

            cy, cx = np.unravel_index(np.argmax(nimg), img.shape)

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

        img_centers[idx, 0] = cx_final
        img_centers[idx, 1] = cy_final

    # look for outliers and replace by a linear fit to all good ones
    # Ticket #81
    ibad = []
    if nwave > 2:
        c_med = np.median(img_centers, axis=0)
        c_std = np.std(img_centers, axis=0)
        bad   = np.any(np.logical_or(img_centers < (c_med-3*c_std),
                                     img_centers > (c_med+3*c_std)), axis=1)
        ibad  = np.where(bad)[0]
        igood = np.where(np.logical_not(bad))[0]
        if len(ibad) != 0:
            logger.info(f'   ==> found {len(ibad)} outliers. Will replace them with a linear fit.')
            
            idx = np.arange(nwave)

            # x
            lin = np.polyfit(idx[igood], img_centers[igood, 0], 1)
            pol = np.poly1d(lin)
            img_centers[ibad, 0] = pol(idx[ibad])

            # y
            lin = np.polyfit(idx[igood], img_centers[igood, 1], 1)
            pol = np.poly1d(lin)
            img_centers[ibad, 1] = pol(idx[ibad])

    #
    # Generate summary plot
    #
    
    # multi-page PDF to save result
    if save_path is not None:
        pdf = PdfPages(save_path)
    
        for idx, (cwave, img) in enumerate(zip(wave, cube)):
            cx_final = img_centers[idx, 0]
            cy_final = img_centers[idx, 1]
            
            failed = (idx in ibad)
            if failed:
                mcolor = 'r'
                bcolor = 'r'
            else:
                mcolor = 'b'
                bcolor = 'w'
            
            plt.figure('PSF center - imaging', figsize=(8.3, 8))
            plt.clf()

            plt.subplot(111)
            plt.imshow(img/np.nanmax(img), aspect='equal', vmin=1e-6, vmax=1, norm=colors.LogNorm(), 
                       interpolation='nearest', cmap=global_cmap)
            plt.plot([cx_final], [cy_final], marker='D', color=mcolor)
            plt.gca().add_patch(patches.Rectangle((cx-box, cy-box), 2*box, 2*box, ec=bcolor, fc='none'))
            if failed:
                plt.text(cx, cy+box, 'Fit failed', color='r', weight='bold', fontsize='x-small',
                         ha='center', va='bottom')
            plt.title(r'Image #{0} - {1:.0f} nm'.format(idx+1, cwave))

            ext = 1000 / pixel
            plt.xlim(cx_final-ext, cx_final+ext)
            plt.xlabel('x position [pix]')
            plt.ylim(cy_final-ext, cy_final+ext)
            plt.ylabel('y position [pix]')

            plt.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

            pdf.savefig()

        pdf.close()
    
    return img_centers


def star_centers_from_PSF_lss_cube(cube, wave_cube, pixel, save_path=None, logger=_log):
    '''
    Compute star center from PSF LSS spectra (IRDIS LSS)

    Parameters
    ----------
    cube : array_like
        LSS PSF cube

    wave_cube : array_like
        Wavelength values for each field, in nanometers

    pixel : float
        Pixel scale, in mas/pixel

    save_path : str
        Path where to save the fit images. Default is None, which means
        that the plot is not produced

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    psf_centers : array_like
        The star center in each frame and wavelength of the cube
    '''

    # standard parameters
    box = 20

    # prepare plot
    if save_path:
        plt.figure('PSF center - spectro', figsize=(6, 12))
        plt.clf()

    # loop over fiels and wavelengths
    nimg = len(cube)
    psf_centers = np.full((1024, nimg), np.nan)
    for fidx, img in enumerate(cube):
        logger.info('   ==> field {0:2d}/{1:2d}'.format(fidx+1, nimg))

        # remove any NaN
        img = np.nan_to_num(cube[fidx])

        # approximate center
        prof = np.sum(img, axis=0)
        cx_int = np.int(np.argmax(prof))

        # sub-image
        sub = img[:, cx_int-box:cx_int+box]
        xx  = np.arange(2*box)

        # wavelengths for this field
        wave = wave_cube[fidx]

        good = np.where(np.isfinite(wave))[0]
        for widx in good:
            # lambda/D
            loD = wave[widx]*1e-9/8 * 180/np.pi * 3600*1000/pixel

            # current profile
            prof = sub[widx, :]

            # gaussian fit
            imax = np.argmax(prof)

            g_init = models.Gaussian1D(amplitude=prof.max(), mean=imax, stddev=loD) + \
                models.Const1D(amplitude=0)

            fit_g = fitting.LevMarLSQFitter()
            par = fit_g(g_init, xx, prof)

            cx = par[0].mean.value - box + cx_int

            psf_centers[widx, fidx] = cx

        if save_path:
            plt.subplot(1, 2, fidx+1)

            plt.imshow(img/img.max(), aspect='equal', vmin=1e-6, vmax=1, norm=colors.LogNorm(), 
                       interpolation='nearest', cmap=global_cmap)
            plt.plot(psf_centers[:, fidx], range(1024), marker='.', color='dodgerblue', linestyle='none',
                     ms=2, alpha=0.5)

            plt.title(r'Field #{0}'.format(fidx+1))

            ext = 1000 / pixel
            plt.xlim(cx_int-ext, cx_int+ext)
            plt.xlabel('x position [pix]')
            
            plt.ylim(0, 1024)
            if fidx == 0:
                plt.ylabel('y position [pix]')
            else:
                plt.gca().yaxis.set_ticklabels([])

    if save_path:
        plt.subplots_adjust(left=0.15, right=0.98, bottom=0.07, top=0.965, wspace=0.05)
        plt.savefig(save_path)

    return psf_centers


def star_centers_from_waffle_img_cube(cube_cen, wave, waffle_orientation, center_guess, pixel, 
                                      orientation_offset,  high_pass=False, center_offset=(0, 0), 
                                      smooth=0, coro=True, save_path=None, logger=_log):
    '''
    Compute star center from waffle images (IRDIS CI, IRDIS DBI, IFS)

    Parameters
    ----------
    cube_cen : array_like
        IRDIFS waffle cube

    wave : array_like
        Wavelength values, in nanometers

    waffle_orientation : str
        String giving the waffle orientation '+' or 'x'

    center_guess : array
        Estimation of the image center as a function of wavelength. 
        This should be an array of shape nwave*2.

    pixel : float
        Pixel scale, in mas/pixel

    orientation_offset : float
        Field orientation offset, in degrees

    high_pass : bool
        Apply high-pass filter to the image before searching for the
        satelitte spots. Default is False

    smooth : int
        Apply a gaussian smoothing to the images to reduce noise. The
        value is the sigma of the gaussian in pixel.  Default is no
        smoothing

    center_offset : tuple
        Apply an (x,y) offset to the default center position. The offset
        will move the search box of the waffle spots by the amount of
        specified pixels in each direction. Default is no offset

    coro : bool
        Observation was performed with a coronagraph. Default is True

    save_path : str
        Path where to save the fit images. Default is None, which means
        that the plot is not produced

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    spot_centers : array_like
        Centers of each individual spot in each frame of the cube

    spot_dist : array_like
        The 6 possible distances between the different spots

    img_centers : array_like
        The star center in each frame of the cube

    '''

    # standard parameters
    nwave = wave.size
    loD = wave*1e-9/8 * 180/np.pi * 3600*1000/pixel

    # waffle parameters
    freq = 10 * np.sqrt(2) * 0.97
    box = 8
    if waffle_orientation == '+':
        orient = orientation_offset * np.pi / 180
    elif waffle_orientation == 'x':
        orient = orientation_offset * np.pi / 180 + np.pi / 4

    # spot fitting
    xx, yy = np.meshgrid(np.arange(2*box), np.arange(2*box))

    # multi-page PDF to save result
    if save_path is not None:
        pdf = PdfPages(save_path)

    # loop over images
    spot_centers = np.zeros((nwave, 4, 2))
    spot_dist    = np.zeros((nwave, 6))
    img_centers  = np.zeros((nwave, 2))
    for idx, (wave, img) in enumerate(zip(wave, cube_cen)):
        logger.info('   ==> wave {0:2d}/{1:2d} ({2:4.0f} nm)'.format(idx+1, nwave, wave))

        # remove any NaN
        img = np.nan_to_num(img)

        # center guess (+offset)
        cx_int = int(center_guess[idx, 0]) + center_offset[0]
        cy_int = int(center_guess[idx, 1]) + center_offset[1]

        # optional high-pass filter
        if high_pass:
            img = img - ndimage.median_filter(img, 15, mode='mirror')

        # optional smoothing
        if smooth > 0:
            img = ndimage.gaussian_filter(img, smooth)

        # mask for non-coronagraphic observations
        if not coro:
            mask = aperture.disc(cube_cen[0].shape[-1], 5*loD[idx], diameter=False,
                                 center=(cx_int, cy_int), invert=True)
            img *= mask

        # create plot if needed
        if save_path:
            fig = plt.figure('Waffle center - imaging', figsize=(8.3, 8))
            plt.clf()

            if high_pass:
                norm = colors.PowerNorm(gamma=1)
                vmin = -1e-1
                vmax = 1e-1
            else:
                norm = colors.LogNorm()
                vmin = 1e-2
                vmax = 1
            
            col = ['green', 'blue', 'deepskyblue', 'purple']
            ax = fig.add_subplot(111)
            ax.imshow(img/img.max(), aspect='equal', vmin=vmin, vmax=vmax, norm=norm,
                      interpolation='nearest', cmap=global_cmap)
            ax.set_title(r'Image #{0} - {1:.0f} nm'.format(idx+1, wave))
            ax.set_xlabel('x position [pix]')
            ax.set_ylabel('y position [pix]')

        # satelitte spots
        for s in range(4):
            cx = int(cx_int + freq*loD[idx] * np.cos(orient + np.pi/2*s))
            cy = int(cy_int + freq*loD[idx] * np.sin(orient + np.pi/2*s))

            sub = img[cy-box:cy+box, cx-box:cx+box]

            # bounds for fitting: spots slightly outside of the box are allowed
            gbounds = {
                'amplitude': (0.0, None),
                'x_mean': (-2.0, box*2+2),
                'y_mean': (-2.0, box*2+2),
                'x_stddev': (1.0, 20.0),
                'y_stddev': (1.0, 20.0)
            }

            # fit: Gaussian + constant
            imax = np.unravel_index(np.argmax(sub), sub.shape)
            g_init = models.Gaussian2D(amplitude=sub.max(), x_mean=imax[1], y_mean=imax[0],
                                       x_stddev=loD[idx], y_stddev=loD[idx], bounds=gbounds) + \
                                       models.Const2D(amplitude=sub.min())
            fitter = fitting.LevMarLSQFitter()
            par = fitter(g_init, xx, yy, sub)
            fit = par(xx, yy)

            cx_final = cx - box + par[0].x_mean
            cy_final = cy - box + par[0].y_mean

            spot_centers[idx, s, 0] = cx_final
            spot_centers[idx, s, 1] = cy_final

            # plot sattelite spots and fit
            if save_path:
                ax.plot([cx_final], [cy_final], marker='D', color=col[s], zorder=1000)
                ax.add_patch(patches.Rectangle((cx-box, cy-box), 2*box, 2*box, ec='white', fc='none'))

                axs = fig.add_axes((0.17+s*0.2, 0.17, 0.1, 0.1))
                axs.imshow(sub, aspect='equal', vmin=0, vmax=sub.max(), interpolation='nearest', 
                           cmap=global_cmap)
                axs.plot([par[0].x_mean], [par[0].y_mean], marker='D', color=col[s])
                axs.set_xticks([])
                axs.set_yticks([])

                axs = fig.add_axes((0.17+s*0.2, 0.06, 0.1, 0.1))
                axs.imshow(fit, aspect='equal', vmin=0, vmax=sub.max(), interpolation='nearest', 
                           cmap=global_cmap)
                axs.set_xticks([])
                axs.set_yticks([])

        # lines intersection
        intersect = lines_intersect(spot_centers[idx, 0, :], spot_centers[idx, 2, :],
                                    spot_centers[idx, 1, :], spot_centers[idx, 3, :])
        img_centers[idx] = intersect

        # scaling
        spot_dist[idx, 0] = np.sqrt(np.sum((spot_centers[idx, 0, :] - spot_centers[idx, 2, :])**2))
        spot_dist[idx, 1] = np.sqrt(np.sum((spot_centers[idx, 1, :] - spot_centers[idx, 3, :])**2))
        spot_dist[idx, 2] = np.sqrt(np.sum((spot_centers[idx, 0, :] - spot_centers[idx, 1, :])**2))
        spot_dist[idx, 3] = np.sqrt(np.sum((spot_centers[idx, 0, :] - spot_centers[idx, 3, :])**2))
        spot_dist[idx, 4] = np.sqrt(np.sum((spot_centers[idx, 1, :] - spot_centers[idx, 2, :])**2))
        spot_dist[idx, 5] = np.sqrt(np.sum((spot_centers[idx, 2, :] - spot_centers[idx, 3, :])**2))

        # finalize plot
        if save_path:
            ax.plot([spot_centers[idx, 0, 0], spot_centers[idx, 2, 0]],
                    [spot_centers[idx, 0, 1], spot_centers[idx, 2, 1]],
                    color='w', linestyle='dashed', zorder=900)
            ax.plot([spot_centers[idx, 1, 0], spot_centers[idx, 3, 0]],
                    [spot_centers[idx, 1, 1], spot_centers[idx, 3, 1]],
                    color='w', linestyle='dashed', zorder=900)

            ax.plot([intersect[0]], [intersect[1]], marker='+', color='w', ms=15)

            ext = 1000 / pixel
            ax.set_xlim(intersect[0]-ext, intersect[0]+ext)
            ax.set_ylim(intersect[1]-ext, intersect[1]+ext)

            plt.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

            if save_path:
                pdf.savefig()

    if save_path:
        pdf.close()

    return spot_centers, spot_dist, img_centers


def star_centers_from_waffle_lss_cube(cube_cen, cube_sci, wave_cube, center_guess, pixel, high_pass=False,
                                      save_path=None, logger=_log):
    '''
    Compute star center from waffle LSS spectra (IRDIS LSS)

    Parameters
    ----------
    cube_cen : array_like
        LSS waffle cube

    cube_sci : array_like
        Science cube with same DIT as the waffle cube. Can be None

    wave_cube : array_like
        Wavelength values for each field, in nanometers

    center_guess : tupple
        Approximate center of the two fields

    pixel : float
        Pixel scale, in mas/pixel

    high_pass : bool
        Apply high-pass filter to the image before searching for the
        satelitte spots. Default is False

    save_path : str
        Path where to save the fit images. Default is None, which means
        that the plot is not produced

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    spot_centers : array_like
        Centers of spots in each frame and each wavelength of the cube

    spot_dist : array_like
        Distance between the spots in each frame and wavelength of the cube

    img_centers : array_like
        The star center in each frame and wavelength of the cube
    '''

    # standard parameters
    box = 120

    # loop over fiels and wavelengths
    nimg = len(cube_cen)

    # prepare plot
    if save_path:
        plt.figure('Waffle centering - spectro', figsize=(6, 12))
        plt.clf()

    # subtract science cube if provided
    if cube_sci is not None:
        logger.info('   ==> subtract science cube')
        cube_cen -= cube_sci

    spot_centers = np.full((1024, 2, 2), np.nan)
    spot_dist    = np.full((1024, nimg), np.nan)
    img_centers  = np.full((1024, nimg), np.nan)
    for fidx, img in enumerate(cube_cen):
        logger.info('   ==> field {0:2d}/{1:2d}'.format(fidx+1, nimg))

        # remove any NaN
        img = np.nan_to_num(cube_cen[fidx])

        if high_pass:
            img = img - ndimage.median_filter(img, 15, mode='mirror')

        # sub-image
        cx_int = center_guess[fidx, 0]
        sub = img[:, cx_int-box:cx_int+box]
        xx  = np.arange(2*box)

        # wavelengths for this field
        wave = wave_cube[fidx]

        good = np.where(np.isfinite(wave))[0]
        for widx in good:
            # lambda/D
            loD = wave[widx]*1e-9/8 * 180/np.pi * 3600*1000/pixel

            # first waffle
            prof = sub[widx] * (xx < box).astype(np.int)
            imax = np.argmax(prof)
            g_init = models.Gaussian1D(amplitude=prof.max(), mean=imax, stddev=loD) + \
                models.Const1D(amplitude=0)
            fit_g = fitting.LevMarLSQFitter()
            par = fit_g(g_init, xx, prof)

            c0 = par[0].mean.value - box + cx_int

            # second waffle
            prof = sub[widx] * (xx > box).astype(np.int)
            imax = np.argmax(prof)
            g_init = models.Gaussian1D(amplitude=prof.max(), mean=imax, stddev=loD) + \
                models.Const1D(amplitude=0)
            fit_g = fitting.LevMarLSQFitter()
            par = fit_g(g_init, xx, prof)

            c1 = par[0].mean.value - box + cx_int

            spot_centers[widx, fidx, 0] = c0
            spot_centers[widx, fidx, 1] = c1

            spot_dist[widx, fidx] = np.abs(c1-c0)

            img_centers[widx, fidx] = (c0 + c1) / 2
            
        if save_path:            
            if high_pass or (cube_sci is not None):
                norm = colors.PowerNorm(gamma=1)
                vmin = -1e-1
                vmax = 1e-1
            else:
                norm = colors.LogNorm()
                vmin = 1e-5
                vmax = 1
            
            plt.subplot(1, 2, fidx+1)
            plt.imshow(img/img.max(), aspect='equal', vmin=vmin, vmax=vmax, interpolation='nearest',
                       cmap=global_cmap, norm=norm)
            plt.plot(spot_centers[:, fidx, 0], range(1024), marker='.', color='dodgerblue', 
                     linestyle='none', ms=2, alpha=1)
            plt.plot(spot_centers[:, fidx, 1], range(1024), marker='.', color='dodgerblue', 
                     linestyle='none', ms=2, alpha=1)
            plt.plot(img_centers[:, fidx], range(1024), marker='.', color='dodgerblue', 
                     linestyle='none', ms=2, alpha=1)

            plt.title(r'Field #{0}'.format(fidx+1))

            ext = 1000 / pixel
            plt.xlim(cx_int-ext, cx_int+ext)
            plt.xlabel('x position [pix]')
            
            plt.ylim(0, 1024)
            if fidx == 0:
                plt.ylabel('y position [pix]')
            else:
                plt.gca().yaxis.set_ticklabels([])

    if save_path:
        plt.subplots_adjust(left=0.15, right=0.98, bottom=0.07, top=0.965, wspace=0.05)
        plt.savefig(save_path)

    return spot_centers, spot_dist, img_centers
