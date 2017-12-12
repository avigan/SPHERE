import os
import glob
import pandas as pd
import subprocess
import numpy as np
import scipy.ndimage as ndimage
import scipy.interpolate as interp
import scipy.optimize as optim
import shutil
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import configparser

from astropy.io import fits
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages

import vltpf.utils.imutils as imutils
import vltpf.utils.aperture as aperture
import vltpf.transmission as transmission
import vltpf.ReductionPath as ReductionPath
import vltpf.toolbox as toolbox


def compute_detector_flat(raw_flat_files, bpm_files=[], mask_vignetting=True):
    '''
    Compute a master detector flat and associated bad pixel map

    Parameters
    ----------
    raw_flat_files : list
        List of 2 raw flat files

    bpm_files : list
        List of bad pixel map files

    mask_vignetting : bool
        Apply a mask on the flats to compensate the optical
        vignetting. The areas of the detector that are vignetted are
        replaced by a value of 1 in the flats. Default is True
    
    Returns
    -------
    flat : array
        Master detector flat

    bpm : array
        Bad pixel map from flat

    '''

    # read bad pixel maps
    bpm_in = toolbox.compute_bad_pixel_map(bpm_files, dtype=np.uint8)

    # read data
    ff0, hdr0 = fits.getdata(raw_flat_files[0], header=True)
    ff1, hdr1 = fits.getdata(raw_flat_files[1], header=True)

    # flatten if needed
    if ff0.ndim == 3:
        ff0 = np.median(ff0, axis=0)

    if ff1.ndim == 3:
        ff1 = np.median(ff1, axis=0)

    # create master flat
    DIT0 = hdr0['HIERARCH ESO DET SEQ1 DIT']
    DIT1 = hdr1['HIERARCH ESO DET SEQ1 DIT']
    
    if DIT0 > DIT1:
        flat = ff0 - ff1
    else:
        flat = ff1 - ff0

    # bad pixels correction    
    flat = imutils.fix_badpix(flat, bpm_in, npix=12, weight=True)

    # flat = imutils.fix_badpix_vip(flat, bpm_in, box=5)
    flat = imutils.sigma_filter(flat, box=5, nsigma=3, iterate=True)
    flat = imutils.sigma_filter(flat, box=7, nsigma=3, iterate=True)

    # normalized flat
    flat = flat / np.median(flat)

    # additional rounad of bad pixels correction
    bpm = (flat <= 0.9) | (flat >= 1.1)
    bpm = bpm.astype(np.uint8)
    flat = imutils.fix_badpix(flat, bpm, npix=12, weight=True)
    # flat = imutils.fix_badpix_vip(flat, bpm_in, box=5)

    # final products
    flat = flat / np.median(flat)
    
    bpm = (flat <= 0.9) | (flat >= 1.1)
    bpm = bpm.astype(np.uint8)

    # apply IFU mask to avoid "edge effects" in the final images,
    # where the the lenslets are vignetted
    if mask_vignetting:
        package_directory = os.path.dirname(os.path.abspath(__file__))
        ifu_mask = fits.getdata(os.path.join(package_directory, 'data', 'ifu_mask.fits'))
        flat[ifu_mask == 0] = 1
    
    return flat, bpm


def sph_ifs_correct_spectral_xtalk(img):
    '''
    Corrects a IFS frame from the spectral crosstalk

    This routines corrects for the SPHERE/IFS spectral crosstalk at
    small scales and (optionally) at large scales. This correction is
    necessary to correct the signal that is "leaking" between
    lenslets. See Antichi et al. (2009ApJ...695.1042A) for a
    theoretical description of the IFS crosstalk. Some informations
    regarding its correction are provided in Vigan et al. (2015), but
    this procedure still lacks a rigorous description and performance
    analysis.

    Since the correction of the crosstalk involves a convolution by a
    kernel of size 41x41, the values at the edges of the frame depend
    on how you choose to apply the convolution. Current implementation
    is EDGE_TRUNCATE. In other parts of the image (i.e. far from the
    edges), the result is identical to original routine by Dino
    Mesa. Note that in the original routine, the convolution that was
    coded did not treat the edges in a clean way defined
    mathematically. The scipy.ndimage.convolve() function offers
    different possibilities for the edges that are all documented.
    
    Parameters
    ----------
    img : array_like
        Input IFS science frame

    Returns
    -------
    img_corr : array_like
        Science frame corrected from the spectral crosstalk

    '''
    
    # definition of the dimension of the matrix
    sepmax = 20
    dim    = sepmax*2+1
    bfac   = 0.727986/1.8

    # defines a matrix to be used around each pixel
    # (the value of the matrix is lower for greater
    # distances form the center.
    x, y = np.meshgrid(np.arange(dim)-sepmax, np.arange(dim)-sepmax)
    rdist  = np.sqrt(x**2 + y**2)
    kernel = 1 / (1+rdist**3 / bfac**3)
    kernel[(np.abs(x) <= 1) & (np.abs(y) <= 1)] = 0

    # convolution and subtraction
    conv = ndimage.convolve(img, kernel, mode='reflect')
    img_corr = img - conv

    return img_corr


def sph_ifs_fix_badpix(img, bpm):
    '''
    Clean the bad pixels in an IFU image

    Extremely effective routine to remove bad pixels. It goes through
    all bad pixels and fit a line beween the first good pixels
    encountered along the same column as the bad pixel, i.e. along the
    spectral axis of each micro-spectrum. Works very well because as
    zeroth-order the spectrum is very smooth and can be approximated
    by a line over one (or a few) bad pixels.

    Parameters
    ----------
    img : array_like
        The image to be cleaned

    bpm : array_like
        Bad pixel map

    Returns
    -------
    img_clean : array_like
        The cleaned image
    '''

    # copy the original image
    img_clean = img.copy()

    # extension over which the good pixels will be looked for along
    # the spectral direction starting from the bad pixel
    ext = 10

    # remove edges in bad pixel map
    bpm[:ext+1, :] = 0
    bpm[:, :ext+1] = 0
    bpm[-ext-1:, :] = 0
    bpm[:, -ext-1:] = 0

    # use NaN for identifying bad pixels directly in the image 
    img_clean[bpm == 1] = np.nan

    # static indices for searching good pixels and for the linear fit
    idx = np.arange(2*ext+1)
    idx_lh = np.arange(ext)+1

    # loop over bad pixels
    badpix = np.where(bpm == 1)
    for y, x in zip(badpix[0], badpix[1]):
        # extract sub-region along the spectral direction
        sub = img_clean[y-ext:y+ext+1, x]

        # sub-regions "above" and "below" the bad pixel
        sub_low = np.flip(img_clean[y-ext//2:y, x], axis=0)
        sub_hig = img_clean[y+1:y+1+ext//2, x]

        # if any of the two is completely bad: skip
        # occurs only in the vignetted areas
        if np.all(np.isnan(sub_low)) or np.all(np.isnan(sub_hig)):
            continue

        # indices of the first good pixels "above" and "below" the bad pixel
        imin_low = idx_lh[~np.isnan(sub_low)].min()
        imin_hig = idx_lh[~np.isnan(sub_hig)].min()

        # linear fit
        xl = idx[ext-imin_low]
        yl = sub[ext-imin_low]

        xh = idx[ext+imin_hig]
        yh = sub[ext+imin_hig]

        a = (yh - yl) / (xh - xl)
        b = yh - a*xh

        fit = a*idx + b

        # replace bad pixel with the fit
        img_clean[y-imin_low+1:y+imin_hig, x] = fit[ext-imin_low+1:ext+imin_hig]

    # put back original value in regions that could not be corrected
    mask = np.isnan(img_clean)
    img_clean[mask] = img[mask]

    return img_clean    


def wavelength_optimisation(wave_ref, wave_scale, wave_lasers, peak_position_lasers):
    '''
    Wavelength optimisation method, to be used in a minimization routine.

    See Vigan et al. (2015, MNRAS, 454, 129) for details of the
    wavelength recalibration:

    https://ui.adsabs.harvard.edu/#abs/2015MNRAS.454..129V/abstract

    Parameters
    ----------
    wave_ref : float
        Reference wavelength (fitted parameter)

    wave_scale : array_like
        Wavelength scaling values
    
    wave_lasers : array_like
        Real wavelength of the calibration lasers; in nanometers.

    peak_position_lasers : array_like
        Position of the peaks of the lasers, in "spectral channel" number

    Returns
    -------
    Difference between the real wavelengths of the lasers and the
    recalibrated value
    '''

    nwave = wave_scale.size
    idx  = np.arange(nwave, dtype=np.float)
    wave = np.full(nwave, wave_ref) * wave_scale
    intrp_func = interp.interp1d(idx, wave, kind='linear')
    wave_peaks = intrp_func(peak_position_lasers)

    diff = wave_peaks - wave_lasers

    return np.max(np.abs(diff))
    

def fit_peak(x, y, display=False):
    '''
    Fit a Gaussian (with linear trend)

    Parameters
    ----------
    x : array_like
        x values

    y : array_like
        y values

    display : bool
        Display the result of the fit
    
    Returns
    -------
    par    
        Fit parameters: Gaussian amplitude, Gaussian mean, Gaussian
        stddev, line slope, line intercept
    '''
    
    # fit: Gaussian + constant
    g_init = models.Gaussian1D(amplitude=y.max(), mean=x[np.argmax(y)]) + models.Linear1D(slope=0, intercept=0)
    fitter = fitting.LevMarLSQFitter()
    fit = fitter(g_init, x, y)

    if display:
        plt.clf()
        plt.plot(x, y, color='k')
        plt.plot(x, fit(x), color='r')
        plt.tight_layout()
    
    return fit.parameters


class Reduction(object):
    '''
    SPHERE/IFS reduction object
    '''

    ##################################################
    # Class variables
    ##################################################

    # specify for each recipe which other recipes need to have been executed before
    recipe_requirements = {
        'sort_files': [],
        'sort_frames': ['sort_files'],
        'check_files_association': ['sort_files'],
        'sph_ifs_cal_dark': ['sort_files'],
        'sph_ifs_cal_detector_flat': ['sort_files', 'sph_ifs_cal_dark'],
        'sph_ifs_cal_specpos': ['sort_files', 'sph_ifs_cal_dark'],
        'sph_ifs_cal_wave': ['sort_files', 'sph_ifs_cal_dark', 'sph_ifs_cal_specpos'],
        'sph_ifs_cal_ifu_flat': ['sort_files', 'sph_ifs_cal_dark', 'sph_ifs_cal_detector_flat',
                                 'sph_ifs_cal_specpos', 'sph_ifs_cal_wave'],
        'sph_ifs_preprocess_science': ['sort_files', 'sort_frames', 'sph_ifs_cal_dark',
                                       'sph_ifs_cal_detector_flat'],
        'sph_ifs_preprocess_wave': ['sort_files', 'sph_ifs_cal_dark', 'sph_ifs_cal_detector_flat'],
        'sph_ifs_science_cubes': ['sort_files', 'sph_ifs_cal_dark', 'sph_ifs_cal_detector_flat',
                                  'sph_ifs_cal_specpos', 'sph_ifs_cal_wave',
                                  'sph_ifs_preprocess_science', 'sph_ifs_preprocess_wave'],
        'sph_ifs_wavelength_recalibration': ['sort_files', 'sort_frames', 'sph_ifs_preprocess_wave',
                                             'sph_ifs_science_cubes'],    
        'sph_ifs_star_center': ['sort_files', 'sort_frames', 'sph_ifs_science_cubes'],
        'sph_ifs_combine_data': ['sort_files', 'sort_frames', 'sph_ifs_science_cubes',
                                 'sph_ifs_wavelength_recalibration', 'sph_ifs_star_center']
    }
    
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

        # expand path
        path = os.path.expanduser(os.path.join(path, ''))
        
        # zeroth-order reduction validation
        raw = os.path.join(path, 'raw')
        if not os.path.exists(raw):
            raise ValueError('No raw/ subdirectory. {0} is not a valid reduction path!'.format(path))
        
        # init path and name
        self._path = ReductionPath.Path(path)
        self._instrument = 'IFS'
        
        # configuration
        package_directory = os.path.dirname(os.path.abspath(__file__))
        configfile = os.path.join(package_directory, 'instruments', self._instrument+'.ini')
        config = configparser.ConfigParser()
        try:
            config.read(configfile)

            # instrument
            self._pixel = float(config.get('instrument', 'pixel'))
            self._nwave = int(config.get('instrument', 'nwave'))
            self._wave_cal_lasers = [float(w) for w in config.get('calibration', 'wave_cal_lasers').split(',')]

            # reduction
            self._config = dict(config.items('reduction'))
            for key, value in self._config.items():
                if (value == 'True'):
                    self._config[key] = True
                elif (value == 'False'):
                    self._config[key] = False
                else:
                    try:
                        value = int(value)
                        self._config[key] = value
                    except ValueError:
                        pass
        except configparser.Error as e:
            raise ValueError('Error reading configuration file for instrument {0}: {1}'.format(self._instrument, e.message))
        
        # execution of recipes
        self._recipe_execution = {
            'sort_files': False,
            'sort_frames': False,
            'check_files_association': False,
            'sph_ifs_cal_dark': False,
            'sph_ifs_cal_detector_flat': False,
            'sph_ifs_cal_specpos': False,
            'sph_ifs_cal_wave': False,
            'sph_ifs_cal_ifu_flat': False,
            'sph_ifs_preprocess_science': False,
            'sph_ifs_preprocess_wave': False,
            'sph_ifs_science_cubes': False,
            'sph_ifs_wavelength_recalibration': False,
            'sph_ifs_star_center': False,
            'sph_ifs_combine_data': False
        }
        
        # reload any existing data frames
        self.read_info()
    
    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return '<Reduction, instrument={0}, path={1}>'.format(self._instrument, self._path)
    
    ##################################################
    # Properties
    ##################################################
    
    @property
    def instrument(self):
        return self._instrument

    @property
    def pixel(self):
        return self._pixel
    
    @property
    def nwave(self):
        return self._nwave
    
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

    @property
    def config(self):
        return self._config    

    ##################################################
    # Generic class methods
    ##################################################

    def show_config(self):
        '''
        Shows the reduction configuration
        '''

        # dictionary
        dico = self._config

        # silent parameter
        print('{0:<30s}{1}'.format('Parameter', 'Value'))
        print('-'*35)
        key = 'silent'
        print('{0:<30s}{1}'.format(key, dico[key]))

        # pre-processing
        print('-'*35)
        keys = [key for key in dico if key.startswith('preproc')]
        for key in keys:
            print('{0:<30s}{1}'.format(key, dico[key]))

        # centring
        print('-'*35)
        keys = [key for key in dico if key.startswith('center')]
        for key in keys:
            print('{0:<30s}{1}'.format(key, dico[key]))
        
        # combining
        print('-'*35)
        keys = [key for key in dico if key.startswith('combine')]
        for key in keys:
            print('{0:<30s}{1}'.format(key, dico[key]))

        # clean
        print('-'*35)
        keys = [key for key in dico if key.startswith('clean')]
        for key in keys:
            print('{0:<30s}{1}'.format(key, dico[key]))
        print('-'*35)
            
        print()

    
    def init_reduction(self):
        '''
        Sort files and frames, perform sanity check
        '''

        # make sure we have sub-directories
        self._path.create_subdirectories()
        
        self.sort_files()
        self.sort_frames()
        self.check_files_association()
        
    
    def create_static_calibrations(self):
        '''
        Create static calibrations, mainly with esorex
        '''
        
        config = self._config
        
        self.sph_ifs_cal_dark(silent=config['silent'])
        self.sph_ifs_cal_detector_flat(silent=config['silent'])
        self.sph_ifs_cal_specpos(silent=config['silent'])
        self.sph_ifs_cal_wave(silent=config['silent'])
        self.sph_ifs_cal_ifu_flat(silent=config['silent'])
        

    def preprocess_science(self):
        '''
        Collapse and correct raw IFU images
        '''

        config = self._config
        
        self.sph_ifs_preprocess_science(subtract_background=config['preproc_subtract_background'],
                                        fix_badpix=config['preproc_fix_badpix'],
                                        correct_xtalk=config['preproc_fix_badpix'],
                                        collapse_science=config['preproc_collapse_science'],
                                        collapse_type=config['preproc_collapse_type'],
                                        coadd_value=config['preproc_coadd_value'],
                                        collapse_psf=config['preproc_collapse_psf'],
                                        collapse_center=config['preproc_collapse_center'])
        self.sph_ifs_preprocess_wave()


    def process_science(self):
        '''
        Generate (x,y,lambda) cubes, recalibrate wavelength, perform star
        center and combine cubes into final (x,y,time,lambda) cubes
        '''

        config = self._config
        
        self.sph_ifs_science_cubes(silent=config['silent'])
        self.sph_ifs_wavelength_recalibration(high_pass=config['center_high_pass'],
                                              offset=config['center_offset'],
                                              display=config['center_display'],
                                              save=config['center_save'])
        self.sph_ifs_star_center(high_pass=config['center_high_pass'],
                                 offset=config['center_offset'],
                                 display=config['center_display'],
                                 save=config['center_save'])
        self.sph_ifs_combine_data(cpix=config['combine_cpix'],
                                  psf_dim=config['combine_psf_dim'],
                                  science_dim=config['combine_science_dim'],
                                  correct_anamorphism=config['combine_correct_anamorphism'],
                                  nocenter=config['combine_nocenter'],
                                  shift_method=config['combine_shift_method'],
                                  save_scaled=config['combine_save_scaled'])

    
    def clean(self):
        '''
        Clean the reduction directory
        '''
        
        config = self._config
        
        if config['clean']:
            self.sph_ifs_clean(delete_raw=config['clean_delete_raw'],
                               delete_products=config['clean_delete_products'])
        
        
    def full_reduction(self):
        '''
        Performs a full reduction of a data set, from the static
        calibrations to the final (x,y,time,lambda) cubes
        '''
        
        self.init_reduction()
        self.create_static_calibrations()
        self.preprocess_science()
        self.process_science()
        self.clean()

    ##################################################
    # SPHERE/IFS methods
    ##################################################
    
    def read_info(self):
        '''
        Read the files, calibs and frames information from disk

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
            if np.any(files_info['PRO CATG'] == 'IFS_MASTER_DARK'):
                self._recipe_execution['sph_ifs_cal_dark'] = True
            if np.any(files_info['PRO CATG'] == 'IFS_MASTER_DFF'):
                self._recipe_execution['sph_ifs_cal_detector_flat'] = True
            if np.any(files_info['PRO CATG'] == 'IFS_SPECPOS'):
                self._recipe_execution['sph_ifs_cal_specpos'] = True
            if np.any(files_info['PRO CATG'] == 'IFS_WAVECALIB'):
                self._recipe_execution['sph_ifs_cal_wave'] = True
            if np.any(files_info['PRO CATG'] == 'IFS_IFU_FLAT_FIELD'):
                self._recipe_execution['sph_ifs_cal_ifu_flat'] = True
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
        if frames_info is not None:
            wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'WAVE,LAMP')]
            self._recipe_execution['sph_ifs_preprocess_wave'] \
                = os.path.exists(os.path.join(path.preproc, wave_file.index[0]+'_preproc.fits'))

            self._recipe_execution['sph_ifs_wavelength_recalibration'] \
                = os.path.exists(os.path.join(path.products, 'wavelength.fits'))

        if frames_info_preproc is not None:
            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ifs_preprocess_science'] = done
            
            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc_?????'.format(file, idx)
                file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ifs_science_cubes'] = done

            done = True
            files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
                                        (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc_centers'.format(file, idx)
                file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ifs_star_center'] = done

        
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
            raise ValueError('No raw FITS files in reduction path')

        print(' * found {0} FITS files in {1}'.format(len(files), path.raw))

        # read list of keywords
        package_directory = os.path.dirname(os.path.abspath(__file__))
        keywords = []
        file = open(os.path.join(package_directory, 'instruments', 'keywords.dat'), 'r')
        for line in file:
            line = line.strip()
            if line:
                if line[0] != '#':
                    keywords.append(line)
        file.close()

        # short keywords
        keywords_short = keywords.copy()
        for idx in range(len(keywords_short)):
            key = keywords_short[idx]
            if key.find('HIERARCH ESO ') != -1:
                keywords_short[idx] = key[13:]

        # files table
        files_info = pd.DataFrame(index=pd.Index(files, name='FILE'), columns=keywords_short, dtype='float')

        for f in files:
            hdu = fits.open(os.path.join(path.raw, f+'.fits'))
            hdr = hdu[0].header

            for k, sk in zip(keywords, keywords_short):
                files_info.loc[f, sk] = hdr.get(k)

            hdu.close()

        # drop files that are not handled, based on DPR keywords
        files_info.dropna(subset=['DPR TYPE'], inplace=True)
        files_info = files_info[(files_info['DPR CATG'] != 'ACQUISITION') & (files_info['DPR TYPE'] != 'OBJECT,AO')]

        # check instruments
        instru = files_info['SEQ ARM'].unique()
        if len(instru) != 1:
            raise ValueError('Sequence is mixing different instruments: {0}'.format(instru))
        
        # processed column
        files_info.insert(len(files_info.columns), 'PROCESSED', False)
        files_info.insert(len(files_info.columns), 'PRO CATG', ' ')

        # convert times
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)

        # sort by acquisition time
        files_info.sort_values(by='DATE-OBS', inplace=True)
        
        # save files_info
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))    
        self._files_info = files_info

        # update recipe execution
        self._recipe_execution['sort_files'] = True

        
    def sort_frames(self):
        '''
        Extract the frames information from the science files and save
        result in a data frame

        calibs : dataframe
            A data frame with the information on all frames
        '''

        print('Extracting frames information')

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sort_frames', self.recipe_requirements)
        
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
        toolbox.compute_times(frames_info)

        # compute angles (ra, dec, parang)
        toolbox.compute_angles(frames_info)

        # save
        frames_info.to_csv(os.path.join(path.preproc, 'frames.csv'))
        self._frames_info = frames_info

        # update recipe execution
        self._recipe_execution['sort_frames'] = True

        #
        # print some info
        #
        cinfo = frames_info[frames_info['DPR TYPE'] == 'OBJECT']
        if len(cinfo) == 0:
            cinfo = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        
        ra_drot   = cinfo['INS4 DROT2 RA'][0]
        ra_drot_h = np.floor(ra_drot/1e4)
        ra_drot_m = np.floor((ra_drot - ra_drot_h*1e4)/1e2)
        ra_drot_s = ra_drot - ra_drot_h*1e4 - ra_drot_m*1e2
        RA = '{:02.0f}:{:02.0f}:{:02.3f}'.format(ra_drot_h, ra_drot_m, ra_drot_s)
        
        dec_drot  = cinfo['INS4 DROT2 DEC'][0]
        sign = np.sign(dec_drot)
        udec_drot  = np.abs(dec_drot)
        dec_drot_d = np.floor(udec_drot/1e4)
        dec_drot_m = np.floor((udec_drot - dec_drot_d*1e4)/1e2)
        dec_drot_s = udec_drot - dec_drot_d*1e4 - dec_drot_m*1e2
        dec_drot_d *= sign
        DEC = '{:02.0f}:{:02.0f}:{:02.2f}'.format(dec_drot_d, dec_drot_m, dec_drot_s)

        pa_start = cinfo['PARANG'][0]
        pa_end   = cinfo['PARANG'][-1]

        date = str(cinfo['DATE'][0])[0:10]
        
        print(' * Object:      {0}'.format(cinfo['OBJECT'][0]))
        print(' * RA / DEC:    {0} / {1}'.format(RA, DEC))
        print(' * Date:        {0}'.format(date))
        print(' * Instrument:  {0}'.format(cinfo['SEQ ARM'][0]))
        print(' * Derotator:   {0}'.format(cinfo['INS4 DROT2 MODE'][0]))
        print(' * Coronagraph: {0}'.format(cinfo['INS COMB ICOR'][0]))
        print(' * Mode:        {0}'.format(cinfo['INS1 MODE'][0]))
        print(' * Filter:      {0}'.format(cinfo['INS2 COMB IFS'][0]))  
        print(' * DIT:         {0:.2f} sec'.format(cinfo['DET SEQ1 DIT'][0]))
        print(' * NDIT:        {0:.0f}'.format(cinfo['DET NDIT'][0]))
        print(' * Texp:        {0:.2f} min'.format(cinfo['DET SEQ1 DIT'].sum()/60))
        print(' * PA:          {0:.2f}° ==> {1:.2f}° = {2:.2f}°'.format(pa_start, pa_end, np.abs(pa_end-pa_start)))


    def check_files_association(self):
        '''
        Performs the calibration files association as a sanity check

        Warnings and errors are reported at the end. Execution is
        interupted in case of error.
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'check_files_association', self.recipe_requirements)
        
        print('Performing file association for calibrations')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # instrument arm
        arm = files_info['SEQ ARM'].unique()
        if len(arm) != 1:
            raise ValueError('Sequence is mixing different instruments: {0}'.format(arm))
        
        # IFS obs mode
        modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()
        if len(modes) != 1:
            raise ValueError('Sequence is mixing YJ and YJH observations.')

        mode = modes[0]
        if mode == 'OBS_YJ':
            mode_short = 'YJ'
        elif mode == 'OBS_H':
            mode_short = 'YJH'
        else:
            raise ValueError('Unknown IFS mode {0}'.format(mode))

        # specific data frame for calibrations
        # keep static calibrations and sky backgrounds
        calibs = files_info[(files_info['DPR CATG'] == 'CALIB') |
                            ((files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == 'SKY'))]

        ###############################################
        # static calibrations not dependent on science
        ###############################################
        error_flag = 0
        warning_flag = 0

        # white flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
        if len(cfiles) < 2:
            error_flag += 1
            print(' * Error: there should be 2 flat files for white lamp, found {0}'.format(len(cfiles)))
        elif len(cfiles) > 2:
            warning_flag += 1
            print(' * Warning: there should be 2 flat files for white lamp, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # 1020 nm flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
        if len(cfiles) < 2:
            error_flag += 1
            print(' * Error: there should be 2 flat files for 1020 nm filter, found {0}'.format(len(cfiles)))
        elif len(cfiles) > 2:
            warning_flag += 1
            print(' * Warning: there should be 2 flat files for 1020 nm filter, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # 1230 nm flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
        if len(cfiles) < 2:
            error_flag += 1
            print(' * Error: there should be 2 flat files for 1230 nm filter, found {0}'.format(len(cfiles)))
        elif len(cfiles) > 2:
            warning_flag += 1
            print(' * Warning: there should be 2 flat files for 1230 nm filter, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # 1300 nm flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
        if len(cfiles) < 2:
            error_flag += 1
            print(' * Error: there should be 2 flat files for 1300 nm filter, found {0}'.format(len(cfiles)))
        elif len(cfiles) > 2:
            warning_flag += 1
            print(' * Warning: there should be 2 flat files for 1300 nm filter, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # 1550 nm flat (YJH mode only)
        if mode_short == 'YJH':
            cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
            if len(cfiles) < 2:
                error_flag += 1
                print(' * Error: there should be 2 flat files for 1550 nm filter, found {0}'.format(len(cfiles)))
        elif len(cfiles) > 2:
            warning_flag += 1
            print(' * Warning: there should be 2 flat files for 1550 nm filter, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # spectra position
        cfiles = calibs[(calibs['DPR TYPE'] == 'SPECPOS,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
        if len(cfiles) == 0:
            error_flag += 1
            print(' * Error: there should be 1 spectra position file, found none.')
        elif len(cfiles) > 1:
            warning_flag += 1
            print(' * Warning: there should be 1 spectra position file, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # wavelength
        cfiles = calibs[(calibs['DPR TYPE'] == 'WAVE,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
        if len(cfiles) == 0:
            error_flag += 1
            print(' * Error: there should be 1 wavelength calibration file, found none.')
        elif len(cfiles) > 1:
            warning_flag += 1
            print(' * Warning: there should be 1 wavelength calibration file, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # IFU flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
        if len(cfiles) == 0:
            error_flag += 1
            print(' * Error: there should be 1 IFU flat file, found none')
        elif len(cfiles) > 1:
            warning_flag += 1
            print(' * Warning: there should be 1 IFU flat file, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']            
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[2:].index, inplace=True)

        # calibs dark file
        cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                        (calibs['DET SEQ1 DIT'].round(2) == 1.65)]
        if len(cfiles) == 0:
            error_flag += 1
            print(' * Error: there is no dark/background for the basic calibrations (DIT=1.65 sec). ' +
                  'It is mandatory to include one to obtain the best data reduction. ' +
                  'A single dark/background file is sufficient, and it can easily be downloaded ' +
                  'from the ESO archive')

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
                      'usually provide a cleaner data reduction')

        # error reporting
        print('There are {0} warning(s) and {1} error(s) in the classification of files'.format(warning_flag, error_flag))
        if error_flag:
            raise ValueError('There is {0} errors that should be solved before proceeding'.format(error_flag))

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))
        self._files_info = files_info
    
        
    def sph_ifs_cal_dark(self, silent=True):
        '''
        Create the dark and background calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_cal_dark', self.recipe_requirements)
        
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
        DITs  = calibs['DET SEQ1 DIT'].unique().round(2)

        for ctype in types:
            for DIT in DITs:            
                cfiles = calibs[(calibs['DPR TYPE'] == ctype) & (calibs['DET SEQ1 DIT'].round(2) == DIT)]
                files = cfiles.index

                # skip non-existing combinations
                if len(cfiles) == 0:
                    continue

                print(' * {0} with DIT={1:.2f} sec ({2} files)'.format(ctype, DIT, len(cfiles)))

                # create sof
                sof = os.path.join(path.sof, 'dark_DIT={0:.2f}.sof'.format(DIT))
                file = open(sof, 'w')
                for f in files:
                    file.write('{0}{1}.fits     {2}\n'.format(path.raw, f, 'IFS_DARK_RAW'))
                file.close()

                # products
                if ctype == 'SKY':
                    loc = 'sky'
                else:
                    loc = 'internal'
                dark_file = 'dark_{0}_DIT={1:.2f}'.format(loc, DIT)
                bpm_file  = 'dark_{0}_bpm_DIT={1:.2f}'.format(loc, DIT)

                # esorex parameters    
                args = ['esorex',
                        '--no-checksum=TRUE',
                        '--no-datamd5=TRUE',
                        'sph_ifs_master_dark',
                        '--ifs.master_dark.coll_alg=2',
                        '--ifs.master_dark.sigma_clip=3.0',
                        '--ifs.master_dark.smoothing=5',
                        '--ifs.master_dark.min_acceptable=0.0',
                        '--ifs.master_dark.max_acceptable=2000.0',
                        '--ifs.master_dark.outfilename={0}{1}.fits'.format(path.calib, dark_file),
                        '--ifs.master_dark.badpixfilename={0}{1}.fits'.format(path.calib, bpm_file),
                        sof]

                # check esorex
                if shutil.which('esorex') is None:
                    raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                                    'that the ESO pipeline is properly installed before running VLTPF.')

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
                files_info.loc[dark_file, 'INS2 MODE'] = cfiles['INS2 MODE'][0]
                files_info.loc[dark_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
                files_info.loc[dark_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
                files_info.loc[dark_file, 'PROCESSED'] = True
                files_info.loc[dark_file, 'PRO CATG'] = 'IFS_MASTER_DARK'

                files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
                files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
                files_info.loc[bpm_file, 'INS2 MODE'] = cfiles['INS2 MODE'][0]
                files_info.loc[bpm_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
                files_info.loc[bpm_file, 'PROCESSED'] = True
                files_info.loc[bpm_file, 'PRO CATG']  = 'IFS_STATIC_BADPIXELMAP'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_dark'] = True


    def sph_ifs_cal_detector_flat(self, silent=True):
        '''
        Create the detector flat calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_cal_detector_flat', self.recipe_requirements)
        
        print('Creating flats')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        raw = files_info[np.logical_not(files_info['PROCESSED'])]
        calibs = raw[(files_info['DPR TYPE'] == 'FLAT,LAMP') | (files_info['DPR TECH'] == 'IMAGE')]

        # IFS obs mode
        mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]    
        if mode == 'OBS_YJ':
            mode_short = 'YJ'
        elif mode == 'OBS_H':
            mode_short = 'YJH'
        else:
            raise ValueError('Unknown IFS mode {0}'.format(mode))

        # bpm files
        cfiles = files_info[files_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP'].index
        bpm_files = [os.path.join(path.calib, f+'.fits') for f in cfiles]

        # loop on wavelengths
        waves = [         0,        1020,        1230,        1300,        1550]
        combs = ['CAL_BB_2', 'CAL_NB1_1', 'CAL_NB2_1', 'CAL_NB3_1', 'CAL_NB4_2']
        lamps = [         5,           1,           2,           3,           4]

        for wave, comb, lamp in zip(waves, combs, lamps):
            print(' * flat for wavelength {0} nm (filter {1}, lamp {2})'.format(wave, comb, lamp))

            cfiles = calibs[calibs['INS2 COMB IFS'] == '{0}_{1}'.format(comb, mode_short)]
            files = [os.path.join(path.raw, f+'.fits') for f in cfiles.index]

            if len(files) == 0:
                continue
            elif len(files) != 2:
                raise ValueError('There should be exactly 2 raw flat files. Found {0}.'.format(len(files)))

            # create the flat and bpm
            flat, bpm = compute_detector_flat(files, bpm_files=bpm_files, mask_vignetting=True)

            # products
            if wave == 0:
                wav = 'white'
            else:
                wav = str(int(wave))
            flat_file = 'master_detector_flat_{0}_l{1}'.format(wav, lamp)
            bpm_file  = 'dff_badpixelname_{0}_l{1}'.format(wav, lamp)        

            hdu = fits.open(os.path.join(path.raw, files[0]))
            fits.writeto(os.path.join(path.calib, flat_file+'.fits'), flat, header=hdu[0].header, output_verify='silentfix', overwrite=True)
            fits.writeto(os.path.join(path.calib, bpm_file+'.fits'), bpm, header=hdu[0].header, output_verify='silentfix', overwrite=True)
            hdu.close()

            # store products
            files_info.loc[flat_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            files_info.loc[flat_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            files_info.loc[flat_file, 'INS2 MODE'] = cfiles['INS2 MODE'][0]
            files_info.loc[flat_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
            files_info.loc[flat_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
            files_info.loc[flat_file, 'PROCESSED'] = True
            files_info.loc[flat_file, 'PRO CATG'] = 'IFS_MASTER_DFF'

            files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            files_info.loc[bpm_file, 'INS2 MODE'] = cfiles['INS2 MODE'][0]
            files_info.loc[bpm_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
            files_info.loc[bpm_file, 'PROCESSED'] = True
            files_info.loc[bpm_file, 'PRO CATG']  = 'IFS_STATIC_BADPIXELMAP'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_detector_flat'] = True


    def sph_ifs_cal_specpos(self, silent=True):
        '''
        Create the specpos calibration

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_cal_specpos', self.recipe_requirements)
        
        print('Creating specpos')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        specpos_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'SPECPOS,LAMP')]
        if len(specpos_file) != 1:
            raise ValueError('There should be exactly 1 raw specpos files. Found {0}.'.format(len(specpos_file)))

        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DARK') & 
                                (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == 1.65)]
        if len(dark_file) == 0:
            raise ValueError('There should at least 1 dark file for calibrations. Found none.')

        # IFS obs mode
        mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]    
        if mode == 'OBS_YJ':
            Hmode = 'FALSE'
        elif mode == 'OBS_H':
            Hmode = 'TRUE'
        else:
            raise ValueError('Unknown IFS mode {0}'.format(mode))

        # create sof
        sof = os.path.join(path.sof, 'specpos.sof')
        file = open(sof, 'w')
        file.write('{0}{1}.fits     {2}\n'.format(path.raw, specpos_file.index[0], 'IFS_SPECPOS_RAW'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, dark_file.index[0], 'IFS_MASTER_DARK'))
        file.close()

        # products
        specp_file = 'spectra_positions'
        
        # esorex parameters    
        args = ['esorex',
                '--no-checksum=TRUE',
                '--no-datamd5=TRUE',
                'sph_ifs_spectra_positions',
                '--ifs.spectra_positions.hmode={0}'.format(Hmode),
                '--ifs.spectra_positions.outfilename={0}{1}.fits'.format(path.calib, specp_file),
                sof]

        # check esorex
        if shutil.which('esorex') is None:
            raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                            'that the ESO pipeline is properly installed before running VLTPF.')
        
        # execute esorex
        if silent:
            proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
        else:
            proc = subprocess.run(args, cwd=path.tmp)

        if proc.returncode != 0:
            raise ValueError('esorex process was not successful')

        # store products
        files_info.loc[specp_file, 'DPR CATG'] = specpos_file['DPR CATG'][0]
        files_info.loc[specp_file, 'DPR TYPE'] = specpos_file['DPR TYPE'][0]
        files_info.loc[specp_file, 'INS2 MODE'] = specpos_file['INS2 MODE'][0]
        files_info.loc[specp_file, 'INS2 COMB IFS'] = specpos_file['INS2 COMB IFS'][0]
        files_info.loc[specp_file, 'DET SEQ1 DIT'] = specpos_file['DET SEQ1 DIT'][0]
        files_info.loc[specp_file, 'PROCESSED'] = True
        files_info.loc[specp_file, 'PRO CATG'] = 'IFS_SPECPOS'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_specpos'] = True


    def sph_ifs_cal_wave(self, silent=True):
        '''
        Create the wavelength calibration

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_cal_wave', self.recipe_requirements)
        
        print('Creating wavelength calibration')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'WAVE,LAMP')]
        if len(wave_file) != 1:
            raise ValueError('There should be exactly 1 raw wavelength calibration file. Found {0}.'.format(len(wave_file)))

        specpos_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_SPECPOS')]
        if len(specpos_file) != 1:
            raise ValueError('There should be exactly 1 specpos file. Found {0}.'.format(len(specpos_file)))

        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DARK') & 
                                (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == 1.65)]
        if len(dark_file) == 0:
            raise ValueError('There should at least 1 dark file for calibrations. Found none.')

        # IFS obs mode
        mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]            

        # create sof
        sof = os.path.join(path.sof, 'wave.sof')
        file = open(sof, 'w')
        file.write('{0}{1}.fits     {2}\n'.format(path.raw, wave_file.index[0], 'IFS_WAVECALIB_RAW'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, specpos_file.index[0], 'IFS_SPECPOS'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, dark_file.index[0], 'IFS_MASTER_DARK'))
        file.close()

        # products
        wav_file = 'wave_calib'
        
        # esorex parameters
        if mode == 'OBS_YJ':
            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ifs_wave_calib',
                    '--ifs.wave_calib.number_lines=3',
                    '--ifs.wave_calib.wavelength_line1=0.9877',
                    '--ifs.wave_calib.wavelength_line2=1.1237',
                    '--ifs.wave_calib.wavelength_line3=1.3094',
                    '--ifs.wave_calib.outfilename={0}{1}.fits'.format(path.calib, wav_file),
                    sof]
        elif mode == 'OBS_H':
            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ifs_wave_calib',
                    '--ifs.wave_calib.number_lines=3',
                    '--ifs.wave_calib.wavelength_line1=0.9877',
                    '--ifs.wave_calib.wavelength_line2=1.1237',
                    '--ifs.wave_calib.wavelength_line3=1.3094',
                    '--ifs.wave_calib.wavelength_line4=1.5451',
                    '--ifs.wave_calib.outfilename={0}{1}.fits'.format(path.calib, wav_file),
                    sof]

        # check esorex
        if shutil.which('esorex') is None:
            raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                            'that the ESO pipeline is properly installed before running VLTPF.')
        
        # execute esorex
        if silent:
            proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
        else:
            proc = subprocess.run(args, cwd=path.tmp)

        if proc.returncode != 0:
            raise ValueError('esorex process was not successful')

        # store products
        files_info.loc[wav_file, 'DPR CATG'] = wave_file['DPR CATG'][0]
        files_info.loc[wav_file, 'DPR TYPE'] = wave_file['DPR TYPE'][0]
        files_info.loc[wav_file, 'INS2 MODE'] = wave_file['INS2 MODE'][0]
        files_info.loc[wav_file, 'INS2 COMB IFS'] = wave_file['INS2 COMB IFS'][0]
        files_info.loc[wav_file, 'DET SEQ1 DIT'] = wave_file['DET SEQ1 DIT'][0]
        files_info.loc[wav_file, 'PROCESSED'] = True
        files_info.loc[wav_file, 'PRO CATG'] = 'IFS_WAVECALIB'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_wave'] = True


    def sph_ifs_cal_ifu_flat(self, silent=True):
        '''
        Create the IFU flat calibration

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_cal_ifu_flat', self.recipe_requirements)
        
        print('Creating IFU flat')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # IFS obs mode
        mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]            
        if mode == 'OBS_YJ':
            mode_short = 'YJ'
        elif mode == 'OBS_H':
            mode_short = 'YJH'
        else:
            raise ValueError('Unknown IFS mode {0}'.format(mode))

        # get list of files
        ifu_flat_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'FLAT,LAMP') &
                                   (files_info['DPR TECH'] == 'IFU')]
        if len(ifu_flat_file) != 1:
            raise ValueError('There should be exactly 1 raw IFU flat file. Found {0}.'.format(len(ifu_flat_file)))

        wave_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_WAVECALIB')]
        if len(wave_file) != 1:
            raise ValueError('There should be exactly 1 wavelength calibration file. Found {0}.'.format(len(wave_file)))

        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DARK') &
                                (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == 1.65)]
        if len(dark_file) == 0:
            raise ValueError('There should at least 1 dark file for calibrations. Found none.')

        flat_white_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                      (files_info['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
        if len(flat_white_file) != 1:
            raise ValueError('There should be exactly 1 white flat file. Found {0}.'.format(len(flat_white_file)))

        flat_1020_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
        if len(flat_1020_file) != 1:
            raise ValueError('There should be exactly 1 1020 nm flat file. Found {0}.'.format(len(flat_1020_file)))

        flat_1230_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
        if len(flat_1230_file) != 1:
            raise ValueError('There should be exactly 1 1230 nm flat file. Found {0}.'.format(len(flat_1230_file)))

        flat_1300_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
        if len(flat_1300_file) != 1:
            raise ValueError('There should be exactly 1 1300 nm flat file. Found {0}.'.format(len(flat_1300_file)))

        if mode == 'OBS_H':
            flat_1550_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                         (files_info['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
            if len(flat_1550_file) != 1:
                raise ValueError('There should be exactly 1 1550 nm flat file. Found {0}.'.format(len(flat_1550_file)))

        # create sof
        sof = os.path.join(path.sof, 'ifu_flat.sof')
        file = open(sof, 'w')
        file.write('{0}{1}.fits     {2}\n'.format(path.raw, ifu_flat_file.index[0], 'IFS_FLAT_FIELD_RAW'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, wave_file.index[0], 'IFS_WAVECALIB'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, dark_file.index[0], 'IFS_MASTER_DARK'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_white_file.index[0], 'IFS_MASTER_DFF_SHORT'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_white_file.index[0], 'IFS_MASTER_DFF_LONGBB'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1020_file.index[0], 'IFS_MASTER_DFF_LONG1'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1230_file.index[0], 'IFS_MASTER_DFF_LONG2'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1300_file.index[0], 'IFS_MASTER_DFF_LONG3'))
        if mode == 'OBS_H':
            file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1550_file.index[0], 'IFS_MASTER_DFF_LONG4'))
        file.close()

        # products
        ifu_file = 'ifu_flat'

        # esorex parameters
        args = ['esorex',
                '--no-checksum=TRUE',
                '--no-datamd5=TRUE',
                'sph_ifs_instrument_flat',
                '--ifs.instrument_flat.nofit=TRUE',
                '--ifs.instrument_flat.ifu_filename={0}{1}.fits'.format(path.calib, ifu_file),
                sof]

        # check esorex
        if shutil.which('esorex') is None:
            raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                            'that the ESO pipeline is properly installed before running VLTPF.')
        
        # execute esorex
        if silent:
            proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
        else:
            proc = subprocess.run(args, cwd=path.tmp)

        if proc.returncode != 0:
            raise ValueError('esorex process was not successful')

        # store products
        files_info.loc[ifu_file, 'DPR CATG'] = ifu_flat_file['DPR CATG'][0]
        files_info.loc[ifu_file, 'DPR TYPE'] = ifu_flat_file['DPR TYPE'][0]
        files_info.loc[ifu_file, 'INS2 MODE'] = ifu_flat_file['INS2 MODE'][0]
        files_info.loc[ifu_file, 'INS2 COMB IFS'] = ifu_flat_file['INS2 COMB IFS'][0]
        files_info.loc[ifu_file, 'DET SEQ1 DIT'] = ifu_flat_file['DET SEQ1 DIT'][0]
        files_info.loc[ifu_file, 'PROCESSED'] = True
        files_info.loc[ifu_file, 'PRO CATG'] = 'IFS_IFU_FLAT_FIELD'

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ifs_cal_ifu_flat'] = True


    def sph_ifs_preprocess_science(self,
                                   subtract_background=True, fix_badpix=True, correct_xtalk=True,
                                   collapse_science=False, collapse_type='mean', coadd_value=2,
                                   collapse_psf=True, collapse_center=True):
        '''Pre-processes the science frames.

        This function can perform multiple steps:
          - collapse of the frames according to different schemes
          - subtract the background
          - correct bad pixels
          - correct the spectral crosstalk

        For the science, 2 collapse methods are available: mean or
        coadd. With mean, the full cubes are mean-combined into a single
        frame. With coadd, the frames are coadded following the
        coadd_value. This can result in lost frames if the number of NDIT
        is not a multiple of coadd_value.

        For the PSFs and star center frames, there is either no collapse
        or a mean collapse.

        The pre-processed frames are saved in the preproc
        sub-directory and will be combined after the (x,y,lambda) cube
        will be created with esorex.

        Parameters
        ----------
        subtract_background : bool
            Performs background subtraction. Default is True

        fix_badpix : bool
            Performs correction of bad pixels. Default is True

        correct_xtalk : bool
            Performs spectral crosstalk correction. Default is True

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
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_preprocess_science', self.recipe_requirements)
        
        print('Pre-processing science files')

        # parameters
        path = self._path
        files_info = self._files_info
        frames_info = self._frames_info
        
        # clean before we start
        files = glob.glob(os.path.join(path.preproc, '*_DIT???_preproc.fits'))
        for file in files:
            os.remove(file)

        # bpm
        if fix_badpix:
            bpm_files = files_info[files_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP'].index
            bpm_files = [os.path.join(path.calib, f+'.fits') for f in bpm_files]

            bpm = toolbox.compute_bad_pixel_map(bpm_files)

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
                        dfiles = files_info[(files_info['PRO CATG'] == 'IFS_MASTER_DARK') &
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
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean')
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none')
                    elif (typ == 'OBJECT,FLUX'):
                        if collapse_psf:
                            print('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean')
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none')
                    elif (typ == 'OBJECT'):
                        if collapse_science:                        
                            if collapse_type == 'mean':
                                print('   ==> collapse: mean ({0} -> 1 frame, 0 dropped)'.format(len(img)))
                                img = np.mean(img, axis=0, keepdims=True)

                                frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean')
                            elif collapse_type == 'coadd':
                                if (not isinstance(coadd_value, int)) or (coadd_value <= 1):
                                    raise TypeError('coadd_value must be an integer >1')

                                coadd_value = int(coadd_value)
                                NDIT = len(img)
                                NDIT_new = NDIT // coadd_value
                                dropped = NDIT % coadd_value

                                if coadd_value > NDIT:
                                    raise ValueError('coadd_value ({0}) must be < NDIT ({1})'.format(coadd_value, NDIT))

                                print('   ==> collapse: coadd by {0} ({1} -> {2} frames, {3} dropped)'.format(coadd_value, NDIT, NDIT_new, dropped))

                                # coadd frames
                                nimg = np.empty((NDIT_new, 2048, 2048), dtype=img.dtype)
                                for f in range(NDIT_new):
                                    nimg[f] = np.mean(img[f*coadd_value:(f+1)*coadd_value], axis=0)
                                img = nimg

                                frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'coadd', coadd_value=coadd_value)
                            else:
                                raise ValueError('Unknown collapse type {0}'.format(collapse_type))
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none')

                    # merge collapse collapsed frames_info
                    frames_info_preproc = pd.concat((frames_info_preproc, frames_info_new))

                    # background subtraction
                    if subtract_background:
                        print('   ==> subtract background')
                        for f in range(len(img)):
                            img[f] -= bkg

                    # bad pixels correction
                    if fix_badpix:
                        print('   ==> correct bad pixels')
                        for f in range(len(img)):
                            frame = img[f]

                            frame = imutils.sigma_filter(frame, box=5, nsigma=3, iterate=True)
                            frame = imutils.sigma_filter(frame, box=7, nsigma=3, iterate=True)
                            frame = sph_ifs_fix_badpix(frame, bpm)
                            img[f] = frame

                    # spectral crosstalk correction
                    if correct_xtalk:
                        print('   ==> correct spectral crosstalk')
                        for f in range(len(img)):
                            frame = img[f]
                            frame = sph_ifs_correct_spectral_xtalk(frame)
                            img[f] = frame

                    # check prensence of coordinates
                    # if not, warn user and add fake one: it could be internal source data
                    if hdr.get('HIERARCH ESO TEL TARG ALPHA') is None:
                        print('Warning: no valid coordinates found in header. Adding fake ones to be able to produce (x,y,lambda) datacubes.')
                        
                        hdr['HIERARCH ESO TEL TARG ALPHA'] =  120000.0
                        hdr['HIERARCH ESO TEL TARG DELTA'] = -900000.0

                            
                    # save DITs individually
                    for f in range(len(img)):
                        frame = img[f].squeeze()                    
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
        self._recipe_execution['sph_ifs_preprocess_science'] = True

        
    def sph_ifs_preprocess_wave(self):
        '''
        Pre-processes the wavelength calibration frame for later
        recalibration of the wavelength
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_preprocess_wave', self.recipe_requirements)
        
        # parameters
        path = self._path
        files_info = self._files_info
                
        print('Pre-processing wavelength calibration file')

        # bpm
        bpm_files = files_info[files_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP'].index
        bpm_files = [os.path.join(path.calib, f+'.fits') for f in bpm_files]
        bpm = toolbox.compute_bad_pixel_map(bpm_files)

        # dark
        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DARK') &
                               (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == 1.65)]
        if len(dark_file) == 0:
            raise ValueError('There should at least 1 dark file for calibrations. Found none.')
        bkg = fits.getdata(os.path.join(path.calib, dark_file.index[0]+'.fits'))

        # wavelength calibration
        wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'WAVE,LAMP')]
        if len(wave_file) != 1:
            raise ValueError('There should be exactly 1 raw wavelength calibration file. Found {0}.'.format(len(wave_file)))
        fname = wave_file.index[0]

        # read data
        print(' * {0}'.format(fname))
        print('   ==> read data')
        img, hdr = fits.getdata(os.path.join(path.raw, fname+'.fits'), header=True)

        # collapse
        print('   ==> collapse: mean')
        img = np.mean(img, axis=0, keepdims=False)

        # background subtraction
        print('   ==> subtract background')
        img -= bkg

        # bad pixels correction
        print('   ==> correct bad pixels')
        img = sph_ifs_fix_badpix(img, bpm)

        # spectral crosstalk correction
        print('   ==> correct spectral crosstalk')
        img = sph_ifs_correct_spectral_xtalk(img)

        # add fake coordinates
        hdr['HIERARCH ESO TEL TARG ALPHA'] =  120000.0
        hdr['HIERARCH ESO TEL TARG DELTA'] = -900000.0

        # save
        fits.writeto(os.path.join(path.preproc, fname+'_preproc.fits'), img, hdr,
                     overwrite=True, output_verify='silentfix')

        # update recipe execution
        self._recipe_execution['sph_ifs_preprocess_wave'] = True


    def sph_ifs_science_cubes(self, silent=True):
        '''
        Create the science cubes from the preprocessed frames

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_science_cubes', self.recipe_requirements)
        
        print('Creating the (x,y,lambda) science cubes')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # clean before we start
        files = glob.glob(os.path.join(path.tmp, '*_DIT???_preproc_?????.fits'))
        for file in files:
            os.remove(file)

        # IFS obs mode
        mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]            
        if mode == 'OBS_YJ':
            mode_short = 'YJ'
        elif mode == 'OBS_H':
            mode_short = 'YJH'
        else:
            raise ValueError('Unknown IFS mode {0}'.format(mode))

        # get list of science files
        sci_files = glob.glob(path.preproc+'*_preproc.fits')
        print(' * found {0} pre-processed files'.format(len(sci_files)))

        # get list of calibration files
        bpm_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP') &
                              (files_info['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]    

        ifu_flat_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_IFU_FLAT_FIELD')]
        if len(ifu_flat_file) != 1:
            raise ValueError('There should be exactly 1 IFU flat file. Found {0}.'.format(len(ifu_flat_file)))

        wave_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_WAVECALIB')]
        if len(wave_file) != 1:
            raise ValueError('There should be exactly 1 wavelength calibration file. Found {0}.'.format(len(wave_file)))

        flat_white_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                      (files_info['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
        if len(flat_white_file) != 1:
            raise ValueError('There should be exactly 1 white flat file. Found {0}.'.format(len(flat_white_file)))

        flat_1020_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
        if len(flat_1020_file) != 1:
            raise ValueError('There should be exactly 1 1020 nm flat file. Found {0}.'.format(len(flat_1020_file)))

        flat_1230_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
        if len(flat_1230_file) != 1:
            raise ValueError('There should be exactly 1 1230 nm flat file. Found {0}.'.format(len(flat_1230_file)))

        flat_1300_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (files_info['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
        if len(flat_1300_file) != 1:
            raise ValueError('There should be exactly 1 1300 nm flat file. Found {0}.'.format(len(flat_1300_file)))

        if mode == 'OBS_H':
            flat_1550_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                         (files_info['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
            if len(flat_1550_file) != 1:
                raise ValueError('There should be exactly 1 1550 nm flat file. Found {0}.'.format(len(flat_1550_file)))

        # create sof
        sof = os.path.join(path.sof, 'science.sof')
        file = open(sof, 'w')
        for f in sci_files:
            file.write('{0}     {1}\n'.format(f, 'IFS_SCIENCE_DR_RAW'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, ifu_flat_file.index[0], 'IFS_IFU_FLAT_FIELD'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, wave_file.index[0], 'IFS_WAVECALIB'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_white_file.index[0], 'IFS_MASTER_DFF_SHORT'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_white_file.index[0], 'IFS_MASTER_DFF_LONGBB'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, bpm_file.index[0], 'IFS_STATIC_BADPIXELMAP'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1020_file.index[0], 'IFS_MASTER_DFF_LONG1'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1230_file.index[0], 'IFS_MASTER_DFF_LONG2'))
        file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1300_file.index[0], 'IFS_MASTER_DFF_LONG3'))
        if mode == 'OBS_H':
            file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_1550_file.index[0], 'IFS_MASTER_DFF_LONG4'))
        file.close()

        # esorex parameters
        print(' * starting esorex')
        args = ['esorex',
                '--no-checksum=TRUE',
                '--no-datamd5=TRUE',
                'sph_ifs_science_dr',
                '--ifs.science_dr.use_adi=0',
                '--ifs.science_dr.spec_deconv=FALSE',
                sof]

        # check esorex
        if shutil.which('esorex') is None:
            raise NameError('esorex does not appear to be in your PATH. Please make sure ' +
                            'that the ESO pipeline is properly installed before running VLTPF.')
        
        # execute esorex
        if silent:
            proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
        else:
            proc = subprocess.run(args, cwd=path.tmp)

        if proc.returncode != 0:
            # raise ValueError('esorex process was not successful')
            print('Error: esorex was not successful. Trying to process some of the frames...')

        # post-process
        print(' * post-processing files')
        files = glob.glob(path.tmp+'*_preproc_*.fits')
        for f in files:
            # read and save only primary extension
            data, header = fits.getdata(f, header=True)
            fits.writeto(f, data, header, overwrite=True, output_verify='silentfix')

        # move files to final directory
        for file in files:
            shutil.move(file, os.path.join(path.preproc, os.path.basename(file)))

        # update recipe execution
        self._recipe_execution['sph_ifs_science_cubes'] = True


    def sph_ifs_wavelength_recalibration(self, high_pass=False, offset=(0, 0), display=False, save=True):
        '''Performs a recalibration of the wavelength, is star center frames
        are available

        See Vigan et al. (2015, MNRAS, 454, 129) for details of the
        wavelength recalibration:

        https://ui.adsabs.harvard.edu/#abs/2015MNRAS.454..129V/abstract

        Parameters
        ----------
        high_pass : bool
            Apply high-pass filter to the image before searching for the satelitte spots

        offset : tuple
            Apply an (x,y) offset to the default center position, for the waffle centering.
            Default is no offset
        
        display : bool
            Display the fit of the satelitte spots. Default is False.

        save : bool
            Save the fit of the sattelite spot for quality check. Default is True,
            although it is a bit slow.

        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_wavelength_recalibration', self.recipe_requirements)
        
        print('Recalibrating wavelength')

        # parameters
        path = self._path
        nwave = self._nwave
        files_info = self._files_info
        frames_info = self._frames_info_preproc

        #
        # DRH wavelength
        #
        print(' * extracting calibrated wavelength')

        # get header of any science file
        science_files = frames_info[frames_info['DPR CATG'] == 'SCIENCE'].index[0]
        fname = '{0}_DIT{1:03d}_preproc_'.format(science_files[0], science_files[1])
        files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
        hdr = fits.getheader(files[0])

        wave_min = hdr['HIERARCH ESO DRS IFS MIN LAMBDA']
        wave_max = hdr['HIERARCH ESO DRS IFS MAX LAMBDA']
        wave_drh = np.linspace(wave_min, wave_max, nwave)

        #
        # star center
        #
        print(' * fitting satelitte spots')

        # get first DIT of first OBJECT,CENTER in the sequence
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) == 0:
            print(' ==> no OBJECT,CENTER file in the data set. Wavelength cannot be recalibrated. ' +
                  'The standard wavelength calibrated by the ESO pripeline will be used.')
            fits.writeto(os.path.join(path.products, 'wavelength.fits'), wave_drh, overwrite=True)

            return wave_drh

        ifs_mode = starcen_files['INS2 COMB IFS'].values[0]
        fname = '{0}_DIT{1:03d}_preproc_'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])    

        files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
        cube, hdr = fits.getdata(files[0], header=True)

        # compute centers from waffle spots
        waffle_orientation = hdr['HIERARCH ESO OCS WAFFLE ORIENT']
        if save:
            save_path = os.path.join(path.products, fname+'spots_fitting.pdf')
        else:
            save_path = None
        spot_center, spot_dist, img_center \
            = toolbox.star_centers_from_waffle_cube(cube, wave_drh, 'IFS', waffle_orientation,
                                                    high_pass=high_pass, center_offset=offset,
                                                    display=display, save_path=save_path)

        # final scaling
        wave_scales = spot_dist / np.full((nwave, 6), spot_dist[0])
        wave_scale  = wave_scales.mean(axis=1)

        #
        # wavelength recalibration
        #    
        print(' * recalibration')

        # find wavelength calibration file name
        wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'WAVE,LAMP')].index[0]
        fname = '{0}_preproc_'.format(wave_file)
        file = glob.glob(os.path.join(path.preproc, fname+'*.fits'))

        # read cube and measure mean flux in all channels
        cube, hdr = fits.getdata(file[0], header=True)
        wave_flux = np.zeros(nwave)
        aper = aperture.disc(cube.shape[-1], 100, diameter=True)
        mask = aper != 0
        for w, f in enumerate(cube):
            wave_flux[w] = f[mask].mean()

        # fit
        wave_idx = np.arange(nwave, dtype=np.float)
        peak_position_lasers = []
        if ifs_mode == 'OBS_YJ':
            # peak 1
            sub_idx  = wave_idx[0:11]
            sub_flux = wave_flux[0:11]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # peak 2
            sub_idx  = wave_idx[10:27]
            sub_flux = wave_flux[10:27]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # peak 3
            sub_idx  = wave_idx[26:]
            sub_flux = wave_flux[26:]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # wavelengths
            wave_lasers = self._wave_cal_lasers[0:3]
        elif ifs_mode == 'OBS_H':
            # peak 1
            sub_idx  = wave_idx[0:8]
            sub_flux = wave_flux[0:8]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # peak 2
            sub_idx  = wave_idx[5:17]
            sub_flux = wave_flux[5:17]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # peak 3
            sub_idx  = wave_idx[14:26]
            sub_flux = wave_flux[14:26]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # peak 4
            sub_idx  = wave_idx[25:]
            sub_flux = wave_flux[25:]
            par = fit_peak(sub_idx, sub_flux)
            peak_position_lasers.append(par[1])

            # wavelengths
            wave_lasers = self._wave_cal_lasers[0:4]

        res = optim.minimize(wavelength_optimisation, 0.9, method='Nelder-Mead',
                             args=(wave_scale, wave_lasers, peak_position_lasers))

        wave_final = np.full(nwave, res.x) * wave_scale

        wave_diff = np.abs(wave_final - wave_drh)*1000
        print('   ==> difference with calibrated wavelength: ' +
              'min={0:.1f} nm, max={1:.1f} nm'.format(wave_diff.min(), wave_diff.max()))

        # save
        print(' * saving')
        fits.writeto(os.path.join(path.products, 'wavelength.fits'), wave_final, overwrite=True)

        #
        # summary plot
        #
        fig = plt.figure(1, figsize=(17, 5.5))
        plt.clf()
        ax = fig.add_subplot(131)
        ax.plot(img_center[:, 0], img_center[:, 1], linestyle='none', marker='+')
        ax.set_xlabel('x center [pix]')
        ax.set_ylabel('y center [pix]')
        ax.set_xlim(img_center[:, 0].mean()+np.array([-3, 3]))
        ax.set_ylim(img_center[:, 1].mean()+np.array([-3, 3]))
        ax.set_title('Frames centers')

        ax = fig.add_subplot(132)
        ax.plot(wave_scales, linestyle='dotted')
        ax.plot(wave_scale, color='k', label='Mean')
        ax.set_xlabel('Spectral channel index')
        ax.set_ylabel('Scaling factor')
        ax.set_title('Spectral scaling')
        ax.legend(loc='upper left')

        ax = fig.add_subplot(133)
        ax.plot(wave_drh, wave_flux, linestyle='dotted', color='k', label='Original')
        ax.plot(wave_final, wave_flux, color='r', label='Recalibrated')
        for w in self._wave_cal_lasers:
            ax.axvline(x=w, linestyle='dashed', color='purple')
        ax.set_xlabel(r'Wavelength [$\mu$m]')
        ax.set_ylabel('Flux')
        ax.legend(loc='upper right')
        ax.set_title('Wavelength calibration')
        plt.tight_layout()

        plt.savefig(os.path.join(path.products, 'wavelegnth_recalibration.pdf'))

        # update recipe execution
        self._recipe_execution['sph_ifs_wavelength_recalibration'] = True


    def sph_ifs_star_center(self, high_pass=False, offset=(0, 0), display=False, save=True):
        '''Determines the star center for all frames where a center can be
        determined (OBJECT,CENTER and OBJECT,FLUX)

        Parameters
        ----------
        high_pass : bool
            Apply high-pass filter to the image before searching for the satelitte spots

        offset : tuple
            Apply an (x,y) offset to the default center position, for the waffle centering.
            Default is no offset
        
        display : bool
            Display the fit of the satelitte spots

        save : bool
            Save the fit of the sattelite spot for quality check. Default is True,
            although it is a bit slow.

        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_star_center', self.recipe_requirements)
        
        print('Star centers determination')

        # parameters
        path = self._path
        nwave = self._nwave
        pixel = self._pixel
        frames_info = self._frames_info_preproc
        
        # start with OBJECT,FLUX
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        if len(flux_files) != 0:
            for file, idx in flux_files.index:
                print('  ==> OBJECT,FLUX: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc_'.format(file, idx)    
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube, hdr = fits.getdata(files[0], header=True)

                # wavelength
                wave_min = hdr['HIERARCH ESO DRS IFS MIN LAMBDA']
                wave_max = hdr['HIERARCH ESO DRS IFS MAX LAMBDA']
                wave_drh = np.linspace(wave_min, wave_max, nwave)

                # centers
                if save:
                    save_path = os.path.join(path.products, fname+'PSF_fitting.pdf')
                else:
                    save_path = None
                img_center = toolbox.star_centers_from_PSF_cube(cube, wave_drh, pixel, display=display, save_path=save_path)

                # save
                fits.writeto(os.path.join(path.preproc, fname+'centers.fits'), img_center, overwrite=True)
                print()
        
        # then OBJECT,CENTER
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) != 0:
            for file, idx in starcen_files.index:
                print('  ==> OBJECT,CENTER: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc_'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube, hdr = fits.getdata(files[0], header=True)

                # wavelength
                wave_min = hdr['HIERARCH ESO DRS IFS MIN LAMBDA']
                wave_max = hdr['HIERARCH ESO DRS IFS MAX LAMBDA']
                wave_drh = np.linspace(wave_min, wave_max, nwave)

                # centers
                waffle_orientation = hdr['HIERARCH ESO OCS WAFFLE ORIENT']
                if save:
                    save_path = os.path.join(path.products, fname+'spots_fitting.pdf')
                else:
                    save_path = None
                spot_center, spot_dist, img_center \
                    = toolbox.star_centers_from_waffle_cube(cube, wave_drh, 'IFS', waffle_orientation,
                                                            high_pass=high_pass, center_offset=offset,
                                                            display=display, save_path=save_path)
                
                # save
                fits.writeto(os.path.join(path.preproc, fname+'centers.fits'), img_center, overwrite=True)
                print()

        # update recipe execution
        self._recipe_execution['sph_ifs_star_center'] = True


    def sph_ifs_combine_data(self, cpix=True, psf_dim=80, science_dim=290, correct_anamorphism=True,
                             shift_method='fft', nocenter=False, save_scaled=False):
        '''Combine and save the science data into final cubes

        All types of data are combined independently: PSFs
        (OBJECT,FLUX), star centers (OBJECT,CENTER) and standard
        coronagraphic images (OBJECT). For each type of data, the
        method saves 3 or 4 different files:
        
          - *_cube: the (x,y,time,lambda) cube
        
          - *_parang: the parallactic angle vector
        
          - *_derot: the derotation angles vector. This vector takes
                     into account the parallactic angle and any
                     instrumental pupil offset. This is the values
                     that need to be used for aligning the images with
                     North up and East left.
        
          - *_frames: a csv file with all the information for every
                      frames. There is one line by time step in the
                      data cube.
        
          - *_cube_scaled: the (x,y,time,lambda) cube with images
                           rescaled spectraly. This is useful if you
                           plan to perform spectral differential
                           imaging in your analysis.

        The method also save a frames.csv file with all the
        information extracted the raw files headers.
                
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

        correct_anamorphism : bool
            Correct the optical anamorphism of the instrument. Default
            is True. See user manual for details.

        nocenter : bool        
            Control if images are centered or not before being
            combined. Useful if centering must be done
            afterwards. Default is False. Note that if nocenter is
            True, the save_scaled option is automatically disabled.
        
        shift_method : str
            Method to scaling and shifting the images: fft or interp.
            Default is fft
        
        save_scaled : bool
            Also save the wavelength-rescaled cubes. Makes the process
            much longer. The default is False

        '''
        
        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ifs_combine_data', self.recipe_requirements)
        
        print('Combine science data')

        # parameters
        path = self._path
        nwave = self._nwave
        frames_info = self._frames_info_preproc
        
        # read final wavelength calibration
        fname = os.path.join(path.products, 'wavelength.fits')
        if not os.path.exists(fname):
            raise FileExistsError('Missing wavelength.fits file. ' +
                                  'You must first run the sph_ifs_wavelength_recalibration() method.')    
        wave = fits.getdata(fname)    

        # max images size
        if psf_dim > 290:
            print('Warning: psf_dim cannot be larger than 290 pix. A value of 290 will be used.')
            psf_dim = 290

        if science_dim > 290:
            print('Warning: science_dim cannot be larger than 290 pix. A value of 290 will be used.')
            science_dim = 290

        # centering
        if nocenter:
            print('Warning: images will not be centered. They will just be combined.')
            shift_method = 'roll'
            centers_default = np.full((nwave, 2), 290//2)

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
                print('  ==> file {0}/{1}: {2}, DIT={3}'.format(file_idx+1, len(flux_files), file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc_'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'?????.fits'))
                cube = fits.getdata(files[0])
                centers = fits.getdata(os.path.join(path.preproc, fname+'centers.fits'))

                # mask values outside of IFS FoV
                cube[cube == 0] = np.nan
                
                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave*1000)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                psf_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                psf_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    if nocenter:
                        cx, cy = centers_default[wave_idx, :]
                    else:
                        cx, cy = centers[wave_idx, :]

                    img  = img[:-1, :-1].astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    psf_cube[wave_idx, file_idx] = nimg[:psf_dim, :psf_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = psf_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0059, 1.0011), method='interp')
                        psf_cube[wave_idx, file_idx] = nimg
                    
                    # wavelength-scaled version
                    if save_scaled:
                        nimg = psf_cube[wave_idx, file_idx]
                        psf_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            flux_files.to_csv(os.path.join(path.products, 'psf_frames.csv'))
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
                print('  ==> file {0}/{1}: {2}, DIT={3}'.format(file_idx+1, len(starcen_files), file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc_'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'?????.fits'))
                cube = fits.getdata(files[0])
                centers = fits.getdata(os.path.join(path.preproc, fname+'centers.fits'))

                # mask values outside of IFS FoV
                cube[cube == 0] = np.nan
                
                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave*1000)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                cen_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                cen_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    if nocenter:
                        cx, cy = centers_default[wave_idx, :]
                    else:
                        cx, cy = centers[wave_idx, :]

                    img  = img[:-1, :-1].astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    cen_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = cen_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0059, 1.0011), method='interp')
                        cen_cube[wave_idx, file_idx] = nimg
                    
                    # wavelength-scaled version
                    if save_scaled:
                        nimg = cen_cube[wave_idx, file_idx]
                        cen_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            starcen_files.to_csv(os.path.join(path.products, 'starcenter_frames.csv'))
            fits.writeto(os.path.join(path.products, 'starcenter_cube.fits'), cen_cube, overwrite=True)
            fits.writeto(os.path.join(path.products, 'starcenter_parang.fits'), cen_parang, overwrite=True)
            fits.writeto(os.path.join(path.products, 'starcenter_derot.fits'), cen_derot, overwrite=True)
            if save_scaled:
                fits.writeto(os.path.join(path.products, 'starcenter_cube_scaled.fits'), cen_cube_scaled, overwrite=True)

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
                print('Warning: no OBJECT,CENTER file in the data set. Images cannot be accurately centred. ' +
                      'They will just be combined.')

                centers = centers_default
            else:
                fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                centers = fits.getdata(os.path.join(path.preproc, fname))

            # final center
            if cpix:
                cc = science_dim // 2
            else:
                cc = (science_dim - 1) / 2

            # final arrays
            sci_cube   = np.zeros((nwave, nfiles, science_dim, science_dim))
            sci_parang = np.zeros(nfiles)
            sci_derot  = np.zeros(nfiles)
            if save_scaled:
                sci_cube_scaled = np.zeros((nwave, nfiles, science_dim, science_dim))

            # read and combine files
            for file_idx, (file, idx) in enumerate(object_files.index):
                print('  ==> file {0}/{1}: {2}, DIT={3}'.format(file_idx+1, len(object_files), file, idx))

                # read data
                fname = '{0}_DIT{1:03d}_preproc_'.format(file, idx)
                files = glob.glob(os.path.join(path.preproc, fname+'*.fits'))
                cube = fits.getdata(files[0])

                # mask values outside of IFS FoV
                cube[cube == 0] = np.nan
                
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

                    img  = img[:-1, :-1].astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    sci_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = sci_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0059, 1.0011), method='interp')
                        sci_cube[wave_idx, file_idx] = nimg
                    
                    # wavelength-scaled version
                    if save_scaled:
                        nimg = sci_cube[wave_idx, file_idx]
                        sci_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            object_files.to_csv(os.path.join(path.products, 'science_frames.csv'))
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


    def sph_ifs_clean(self, delete_raw=False, delete_products=False):
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
