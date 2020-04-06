import pandas as pd
import subprocess
import logging
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
import collections

from pathlib import Path
from astropy.io import fits
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages

import sphere
import sphere.utils as utils
import sphere.utils.imutils as imutils
import sphere.utils.aperture as aperture
import sphere.transmission as transmission
import sphere.toolbox as toolbox

_log = logging.getLogger(__name__)


def get_wavelength_calibration(wave_calib, centers, wave_min, wave_max):
    '''
    Return the linear wavelength calibration for each IRDIS field

    Parameters
    ----------
    wave_calib : array
        Wavelength calibration data computed by esorex recipe

    centers : tuple
        Center of each field

    wave_min : float
        Minimal usable wavelength

    wave_max : float
        Maximal usable wavelength

    Returns
    -------
    wave_lin : array
        Array with the linear calibration for each field, as a function
        of pixel coordinate
    '''
    wave_map = np.zeros((2, 1024, 1024))
    wave_map[0] = wave_calib[:, 0:1024]
    wave_map[1] = wave_calib[:, 1024:]
    wave_map[(wave_map < wave_min) | (wave_max < wave_map)] = np.nan

    wave_ext = 10
    wave_lin = np.zeros((2, 1024))

    wave_lin[0] = np.mean(wave_map[0, :, centers[0, 0]-wave_ext:centers[0, 0]+wave_ext], axis=1)
    wave_lin[1] = np.mean(wave_map[1, :, centers[1, 0]-wave_ext:centers[1, 0]+wave_ext], axis=1)

    return wave_lin


class SpectroReduction(object):
    '''
    SPHERE/IRDIS long-slit spectroscopy reduction class. It handles
    both the low and medium resolution modes (LRS, MRS)
    '''

    ##################################################
    # Class variables
    ##################################################

    # specify for each recipe which other recipes need to have been executed before
    recipe_requirements = {
        'sort_files': [],
        'sort_frames': ['sort_files'],
        'check_files_association': ['sort_files'],
        'sph_ird_cal_dark': ['sort_files'],
        'sph_ird_cal_detector_flat': ['sort_files'],
        'sph_ird_cal_wave': ['sort_files', 'sph_ird_cal_detector_flat'],
        'sph_ird_preprocess_science': ['sort_files', 'sort_frames', 'sph_ird_cal_dark',
                                       'sph_ird_cal_detector_flat'],
        'sph_ird_star_center': ['sort_files', 'sort_frames', 'sph_ird_cal_wave'],
        'sph_ird_wavelength_recalibration': ['sort_files', 'sort_frames', 'sph_ird_cal_wave'],
        'sph_ird_combine_data': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science'],
        'sph_ird_clean': []
    }

    ##################################################
    # Constructor
    ##################################################

    def __new__(cls, path, log_level='info', sphere_handler=None):
        '''Custom instantiation for the class and initialization for the
           instances

        The customized instantiation enables to check that the
        provided path is a valid reduction path. If not, None will be
        returned for the reduction being created. Otherwise, an
        instance is created and returned at the end.

        Parameters
        ----------
        path : str
            Path to the directory containing the dataset

        level : {'debug', 'info', 'warning', 'error', 'critical'}
            The log level of the handler

        sphere_handler : log handler
            Higher-level SPHERE.Dataset log handler

        '''

        #
        # make sure we are dealing with a proper reduction directory
        #

        # init path
        path = Path(path).expanduser().resolve()

        # zeroth-order reduction validation
        raw = path / 'raw'
        if not raw.exists():
            _log.error('No raw/ subdirectory. {0} is not a valid reduction path'.format(path))
            return None
        else:
            # it's all good: create instance!
            reduction = super(SpectroReduction, cls).__new__(cls)

        #
        # basic init
        #

        # init path
        reduction._path = utils.ReductionPath(path)

        # instrument and mode
        reduction._instrument = 'IRDIS'
        reduction._mode = 'Unknown'

        #
        # logging
        #
        logger = logging.getLogger(str(path))
        logger.setLevel(log_level.upper())
        if logger.hasHandlers():
            for hdlr in logger.handlers:
                logger.removeHandler(hdlr)

        handler = logging.FileHandler(reduction._path.products / 'reduction.log', mode='w', encoding='utf-8')
        formatter = logging.Formatter('%(asctime)s\t%(levelname)8s\t%(message)s')
        formatter.default_msec_format = '%s.%03d'
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        if sphere_handler:
            logger.addHandler(sphere_handler)
        
        reduction._logger = logger

        reduction._logger.info('Creating IRDIS spectroscopy reduction at path {}'.format(path))

        #
        # configuration
        #
        configfile = f'{Path(sphere.__file__).parent}/instruments/{reduction._instrument}.ini'
        config = configparser.ConfigParser()

        reduction._logger.debug('> read configuration')
        config.read(configfile)

        # instrument
        reduction._pixel = float(config.get('instrument', 'pixel'))
        reduction._nwave = -1

        # calibration
        reduction._wave_cal_lasers = np.array(eval(config.get('calibration', 'wave_cal_lasers')))

        # spectro calibration
        reduction._default_center_lrs = np.array(eval(config.get('calibration-spectro', 'default_center_lrs')))
        reduction._wave_min_lrs = eval(config.get('calibration-spectro', 'wave_min_lrs'))
        reduction._wave_max_lrs = eval(config.get('calibration-spectro', 'wave_max_lrs'))

        reduction._default_center_mrs = np.array(eval(config.get('calibration-spectro', 'default_center_mrs')))
        reduction._wave_min_mrs = eval(config.get('calibration-spectro', 'wave_min_mrs'))
        reduction._wave_max_mrs = eval(config.get('calibration-spectro', 'wave_max_mrs'))

        # reduction parameters
        reduction._config = {}
        for group in ['reduction', 'reduction-spectro']:
            items = dict(config.items(group))
            reduction._config.update(items)
            for key, value in items.items():
                try:
                    val = eval(value)
                except NameError:
                    val = value
                reduction._config[key] = val

        #
        # reduction status
        #
        reduction._status = sphere.INIT
        reduction._recipes_status = collections.OrderedDict()

        # reload any existing data frames
        reduction._read_info()

        #
        # return instance
        #
        return reduction

    ##################################################
    # Representation
    ##################################################

    def __repr__(self):
        return '<SpectroReduction, instrument={}, mode={}, path={}>'.format(self._instrument, self._mode, self._path)

    def __format__(self):
        return self.__repr__()

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
    def recipes_status(self):
        return self._recipes_status

    @property
    def config(self):
        return self._config

    @property
    def mode(self):
        return self._mode

    ##################################################
    # Generic class methods
    ##################################################

    def show_config(self):
        '''
        Shows the reduction configuration
        '''

        # dictionary
        dico = self._config

        # misc parameters
        print()
        print('{0:<30s}{1}'.format('Parameter', 'Value'))
        print('-'*35)
        keys = [key for key in dico if key.startswith('misc')]
        for key in keys:
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

        # wave
        print('-'*35)
        keys = [key for key in dico if key.startswith('wave')]
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

        self._logger.info('====> Init <====')

        # make sure we have sub-directories
        self._path.create_subdirectories()

        self.sort_files()
        self.sort_frames()
        self.check_files_association()


    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        self._logger.info('====> Static calibrations <====')

        config = self._config

        self.sph_ird_cal_dark(silent=config['misc_silent_esorex'])
        self.sph_ird_cal_detector_flat(silent=config['misc_silent_esorex'])
        self.sph_ird_cal_wave(silent=config['misc_silent_esorex'])


    def preprocess_science(self):
        '''
        Clean and collapse images
        '''

        self._logger.info('====> Science pre-processing <====')

        config = self._config

        self.sph_ird_preprocess_science(subtract_background=config['preproc_subtract_background'],
                                        fix_badpix=config['preproc_fix_badpix'],
                                        collapse_science=config['preproc_collapse_science'],
                                        collapse_psf=config['preproc_collapse_psf'],
                                        collapse_center=config['preproc_collapse_center'])


    def process_science(self):
        '''
        Perform star center, combine cubes into final (x,y,time,lambda)
        cubes, correct anamorphism and scale the images
        '''

        self._logger.info('====> Science processing <====')

        config = self._config

        self.sph_ird_star_center(high_pass=config['center_high_pass'],
                                 plot=config['misc_plot'])
        self.sph_ird_wavelength_recalibration(fit_scaling=config['wave_fit_scaling'],
                                              plot=config['misc_plot'])
        self.sph_ird_combine_data(cpix=config['combine_cpix'],
                                  psf_dim=config['combine_psf_dim'],
                                  science_dim=config['combine_science_dim'],
                                  correct_mrs_chromatism=config['combine_correct_mrs_chromatism'],
                                  split_posang=config['combine_split_posang'],
                                  shift_method=config['combine_shift_method'],
                                  manual_center=config['combine_manual_center'],
                                  coarse_centering=config['combine_coarse_centering'])

    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''

        self._logger.info('====> Clean-up <====')

        config = self._config

        if config['clean']:
            self.sph_ird_clean(delete_raw=config['clean_delete_raw'],
                               delete_products=config['clean_delete_products'])


    def full_reduction(self):
        '''
        Performs a full reduction of a data set, from the static
        calibrations to the final (x,y,time,lambda) cubes
        '''

        self._logger.info('====> Full reduction <====')

        self.init_reduction()
        self.create_static_calibrations()
        self.preprocess_science()
        self.process_science()
        self.clean()

    ##################################################
    # Private methods
    ##################################################

    def _read_info(self):
        '''
        Read the files, calibs and frames information from disk

        files_info : dataframe
            The data frame with all the information on files

        frames_info : dataframe
            The data frame with all the information on science frames

        frames_info_preproc : dataframe
            The data frame with all the information on science frames after pre-processing

        This function is not supposed to be called directly by the user.

        '''

        self._logger.info('Read existing reduction information')

        # path
        path = self._path

        # files info
        fname = path.preproc / 'files.csv'
        if fname.exists():
            self._logger.debug('> read files.csv')

            files_info = pd.read_csv(fname, index_col=0)

            # convert times
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

            # update recipe execution
            self._update_recipe_status('sort_files', sphere.SUCCESS)
            if np.any(files_info['PRO CATG'] == 'IRD_MASTER_DARK'):
                self._update_recipe_status('sph_ird_cal_dark', sphere.SUCCESS)
            if np.any(files_info['PRO CATG'] == 'IRD_FLAT_FIELD'):
                self._update_recipe_status('sph_ird_cal_detector_flat', sphere.SUCCESS)
            if np.any(files_info['PRO CATG'] == 'IRD_WAVECALIB'):
                self._update_recipe_status('sph_ird_cal_wave', sphere.SUCCESS)

            # update instrument mode
            self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]
        else:
            files_info = None

        fname = path.preproc / 'frames.csv'
        if fname.exists():
            self._logger.debug('> read frames.csv')

            frames_info = pd.read_csv(fname, index_col=(0, 1))

            # convert times
            frames_info['DATE-OBS'] = pd.to_datetime(frames_info['DATE-OBS'], utc=False)
            frames_info['DATE'] = pd.to_datetime(frames_info['DATE'], utc=False)
            frames_info['DET FRAM UTC'] = pd.to_datetime(frames_info['DET FRAM UTC'], utc=False)
            frames_info['TIME START'] = pd.to_datetime(frames_info['TIME START'], utc=False)
            frames_info['TIME'] = pd.to_datetime(frames_info['TIME'], utc=False)
            frames_info['TIME END'] = pd.to_datetime(frames_info['TIME END'], utc=False)

            # update recipe execution
            self._update_recipe_status('sort_frames', sphere.SUCCESS)
        else:
            frames_info = None

        fname = path.preproc / 'frames_preproc.csv'
        if fname.exists():
            self._logger.debug('> read frames_preproc.csv')

            frames_info_preproc = pd.read_csv(fname, index_col=(0, 1))

            # convert times
            frames_info_preproc['DATE-OBS'] = pd.to_datetime(frames_info_preproc['DATE-OBS'], utc=False)
            frames_info_preproc['DATE'] = pd.to_datetime(frames_info_preproc['DATE'], utc=False)
            frames_info_preproc['DET FRAM UTC'] = pd.to_datetime(frames_info_preproc['DET FRAM UTC'], utc=False)
            frames_info_preproc['TIME START'] = pd.to_datetime(frames_info_preproc['TIME START'], utc=False)
            frames_info_preproc['TIME'] = pd.to_datetime(frames_info_preproc['TIME'], utc=False)
            frames_info_preproc['TIME END'] = pd.to_datetime(frames_info_preproc['TIME END'], utc=False)
        else:
            frames_info_preproc = None

        # save data frames in instance variables
        self._files_info = files_info
        self._frames_info = frames_info
        self._frames_info_preproc = frames_info_preproc

        # additional checks to update recipe execution
        if frames_info_preproc is not None:
            done = (path.preproc / 'wavelength_default.fits').exists()
            if done:
                self._update_recipe_status('sph_ird_cal_wave', sphere.SUCCESS)
            self._logger.debug('> sph_ird_cal_wave status = {}'.format(done))

            done = (path.preproc / 'wavelength_recalibrated.fits').exists()
            if done:
                self._update_recipe_status('sph_ird_wavelength_recalibration', sphere.SUCCESS)
            self._logger.debug('> sph_ird_wavelength_recalibration status = {}'.format(done))

            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                file = list(path.preproc.glob('{}.fits'.format(fname)))
                done = done and (len(file) == 1)
            if done:
                self._update_recipe_status('sph_ird_preprocess_science', sphere.SUCCESS)
            self._logger.debug('> sph_ird_preprocess_science status = {}'.format(done))

            done = True
            files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
                                        (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc_centers'.format(file, idx)
                file = list(path.preproc.glob('{}.fits'.format(fname)))
                done = done and (len(file) == 1)
            if done:
                self._update_recipe_status('sph_ird_star_center', sphere.SUCCESS)
            self._logger.debug('> sph_ird_star_center status = {}'.format(done))

        # reduction status
        self._status = sphere.INCOMPLETE


    def _update_recipe_status(self, recipe, status):
        '''Update execution status for reduction and recipe

        Parameters
        ----------
        recipe : str
            Recipe name

        status : sphere status (int)
            Status of the recipe. Can be either one of sphere.NOTSET,
            sphere.SUCCESS or sphere.ERROR
        '''

        self._logger.debug('> update recipe execution')

        self._recipes_status[recipe] = status
        self._recipes_status.move_to_end(recipe)

    ##################################################
    # SPHERE/IRDIS methods
    ##################################################

    def sort_files(self):
        '''
        Sort all raw files and save result in a data frame

        files_info : dataframe
            Data frame with the information on raw files
        '''

        self._logger.info('Sort raw files')

        # update recipe execution
        self._update_recipe_status('sort_files', sphere.NOTSET)

        # parameters
        path = self._path

        # list files
        files = path.raw.glob('*.fits')
        files = [f.stem for f in files]

        if len(files) == 0:
            self._logger.critical('No raw FITS files in reduction path')
            self._update_recipe_status('sort_files', sphere.ERROR)
            self._status = sphere.FATAL
            return

        self._logger.info(' * found {0} raw FITS files'.format(len(files)))

        # read list of keywords
        self._logger.debug('> read keyword list')
        keywords = []
        file = open(Path(sphere.__file__).parent / 'instruments' / 'keywords.dat', 'r')
        for line in file:
            line = line.strip()
            if line:
                if line[0] != '#':
                    keywords.append(line)
        file.close()

        # short keywords
        self._logger.debug('> translate into short keywords')
        keywords_short = keywords.copy()
        for idx in range(len(keywords_short)):
            key = keywords_short[idx]
            if key.find('HIERARCH ESO ') != -1:
                keywords_short[idx] = key[13:]

        # files table
        self._logger.debug('> create files_info data frame')
        files_info = pd.DataFrame(index=pd.Index(files, name='FILE'), columns=keywords_short, dtype='float')

        self._logger.debug('> read FITS keywords')
        for f in files:
            hdu = fits.open(path.raw / '{}.fits'.format(f))
            hdr = hdu[0].header

            for k, sk in zip(keywords, keywords_short):
                files_info.loc[f, sk] = hdr.get(k)

            hdu.close()

        # drop files that are not handled, based on DPR keywords
        self._logger.debug('> drop unsupported file types')
        files_info.dropna(subset=['DPR TYPE'], inplace=True)
        files_info = files_info[(files_info['DPR CATG'] != 'ACQUISITION') & (files_info['DPR TYPE'] != 'OBJECT,AO')]

        # check instruments
        instru = files_info['SEQ ARM'].unique()
        if len(instru) != 1:
            self._logger.critical('Sequence is mixing different instruments: {0}'.format(instru))
            self._update_recipe_status('sort_files', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # check science files
        sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] != 'SKY')]
        if len(sci_files) == 0:
            self._logger.critical('This dataset contains no science frame. There should be at least one!')
            self._update_recipe_status('sort_frames', sphere.ERROR)
            self._status = sphere.FATAL
            return
        
        # processed column
        files_info.insert(len(files_info.columns), 'PROCESSED', False)
        files_info.insert(len(files_info.columns), 'PRO CATG', ' ')

        # convert times
        self._logger.debug('> convert times')
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

        # update instrument mode
        self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]

        # sort by acquisition time
        files_info.sort_values(by='DATE-OBS', inplace=True)

        # save files_info
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')
        self._files_info = files_info

        # update recipe execution
        self._update_recipe_status('sort_files', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sort_frames(self):
        '''
        Extract the frames information from the science files and save
        result in a data frame

        calibs : dataframe
            A data frame with the information on all frames
        '''

        self._logger.info('Extract frames information')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sort_frames',
                                         self.recipe_requirements, logger=self._logger):
            return

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
        self._logger.debug('> create frames_info data frame')
        frames_info = pd.DataFrame(columns=sci_files.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        frames_info = frames_info.align(files_info, level=0)[1]

        # compute timestamps
        toolbox.compute_times(frames_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(frames_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sort_frames', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # save
        self._logger.debug('> save frames.csv')
        frames_info.to_csv(path.preproc / 'frames.csv')
        self._frames_info = frames_info

        #
        # print some info
        #
        self._logger.debug('> print observation info')
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

        posang   = cinfo['INS4 DROT2 POSANG'].unique()

        date = str(cinfo['DATE'][0])[0:10]

        self._logger.info(' * Object:      {0}'.format(cinfo['OBJECT'][0]))
        self._logger.info(' * RA / DEC:    {0} / {1}'.format(RA, DEC))
        self._logger.info(' * Date:        {0}'.format(date))
        self._logger.info(' * Instrument:  {0}'.format(cinfo['SEQ ARM'][0]))
        self._logger.info(' * Derotator:   {0}'.format(cinfo['INS4 DROT2 MODE'][0]))
        self._logger.info(' * Coronagraph: {0}'.format(cinfo['INS COMB ICOR'][0]))
        self._logger.info(' * Mode:        {0}'.format(cinfo['INS1 MODE'][0]))
        self._logger.info(' * Filter:      {0}'.format(cinfo['INS COMB IFLT'][0]))
        self._logger.info(' * DIT:         {0:.2f} sec'.format(cinfo['DET SEQ1 DIT'][0]))
        self._logger.info(' * NDIT:        {0:.0f}'.format(cinfo['DET NDIT'][0]))
        self._logger.info(' * Texp:        {0:.2f} min'.format(cinfo['DET SEQ1 DIT'].sum()/60))
        self._logger.info(' * PA:          {0:.2f}° ==> {1:.2f}° = {2:.2f}°'.format(pa_start, pa_end, np.abs(pa_end-pa_start)))
        self._logger.info(' * POSANG:      {0}'.format(', '.join(['{:.2f}°'.format(p) for p in posang])))

        # update recipe execution
        self._update_recipe_status('sort_frames', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def check_files_association(self):
        '''
        Performs the calibration files association as a sanity check.

        Warnings and errors are reported at the end. Execution is
        interupted in case of error.
        '''

        self._logger.info('File association for calibrations')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'check_files_association',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        files_info = self._files_info

        # instrument arm
        arm = files_info['SEQ ARM'].unique()
        if len(arm) != 1:
            self._logger.error('Sequence is mixing different instruments: {0}'.format(arm))
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        # IRDIS obs mode and filter combination
        modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'].unique()
        if len(modes) != 1:
            self._logger.eror('Sequence is mixing different types of observations: {0}'.format(modes))
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        filter_combs = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()
        if len(filter_combs) != 1:
            self._logger.error('Sequence is mixing different types of filters combinations: {0}'.format(filter_combs))
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        filter_comb = filter_combs[0]
        if (filter_comb != 'S_LR') and (filter_comb != 'S_MR'):
            self._logger.error('Unknown IRDIS-LSS filter combination/mode {0}'.format(filter_comb))
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        # specific data frame for calibrations
        # keep static calibrations and sky backgrounds
        self._logger.debug('> select calib files')
        calibs = files_info[(files_info['DPR CATG'] == 'CALIB') |
                            ((files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == 'SKY'))]

        ###############################################
        # static calibrations not dependent on DIT
        ###############################################
        error_flag = 0
        warning_flag = 0

        # flat
        self._logger.debug('> check instrument flat requirements')
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS COMB IFLT'] == filter_comb)]
        if len(cfiles) <= 1:
            error_flag += 1
            self._logger.error(' * there should be more than 1 flat in filter combination {0}'.format(filter_comb))

        # wave
        self._logger.debug('> check wavelength calibration requirements')
        cfiles = calibs[(calibs['DPR TYPE'] == 'LAMP,WAVE') & (calibs['INS COMB IFLT'] == filter_comb)]
        if len(cfiles) == 0:
            error_flag += 1
            self._logger.error(' * there should be 1 wavelength calibration file, found none.')
        elif len(cfiles) > 1:
            warning_flag += 1
            self._logger.warning(' * there should be 1 wavelength calibration file, found {0}. Using the closest from science.'.format(len(cfiles)))

            # find the two closest to science files
            sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE')]
            time_sci   = sci_files['DATE-OBS'].min()
            time_flat  = cfiles['DATE-OBS']
            time_delta = np.abs(time_sci - time_flat).argsort()

            # drop the others
            files_info.drop(time_delta[1:].index, inplace=True)

        ##################################################
        # static calibrations that depend on science DIT
        ##################################################

        self._logger.debug('> select science files')
        obj = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'DPR TYPE'].apply(lambda s: s[0:6])
        DITs = files_info.loc[(files_info['DPR CATG'] == 'SCIENCE') & (obj == 'OBJECT'), 'DET SEQ1 DIT'].unique().round(2)

        # handle darks in a slightly different way because there might be several different DITs
        self._logger.debug('> check dark/background requirements')
        for DIT in DITs:
            # instrumental backgrounds
            cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                            (calibs['DET SEQ1 DIT'].round(2) == DIT)]
            if len(cfiles) == 0:
                warning_flag += 1
                self._logger.warning(' * there is no dark/background for science files with DIT={0} sec. It is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive'.format(DIT))

            # sky backgrounds
            cfiles = files_info[(files_info['DPR TYPE'] == 'SKY') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
            if len(cfiles) == 0:
                warning_flag += 1
                self._logger.warning(' * there is no sky background for science files with DIT={0} sec. Using a sky background instead of an internal instrumental background can usually provide a cleaner data reduction'.format(DIT))

        # error reporting
        self._logger.debug('> report status')
        if error_flag:
            self._logger.error('There are {0} warning(s) and {1} error(s) in the classification of files'.format(warning_flag, error_flag))
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return
        else:
            self._logger.warning('There are {0} warning(s) and {1} error(s) in the classification of files'.format(warning_flag, error_flag))

        # save
        files_info.to_csv(path.preproc / 'files.csv')
        self._files_info = files_info

        # update recipe execution
        self._update_recipe_status('check_files_association', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_cal_dark(self, silent=True):
        '''
        Create the dark and background calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        self._logger.info('Darks and backgrounds')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_cal_dark',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        files_info = self._files_info

        # get list of files
        calibs = files_info[np.logical_not(files_info['PROCESSED']) &
                            ((files_info['DPR TYPE'] == 'DARK') |
                             (files_info['DPR TYPE'] == 'DARK,BACKGROUND') |
                             (files_info['DPR TYPE'] == 'SKY'))]

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

                    self._logger.info(' * {0} in filter {1} with DIT={2:.2f} sec ({3} files)'.format(ctype, cfilt, DIT, len(cfiles)))

                    # create sof
                    self._logger.debug('> create sof file')
                    sof = path.sof / 'dark_filt={0}_DIT={1:.2f}.sof'.format(cfilt, DIT)
                    file = open(sof, 'w')
                    for f in files:
                        file.write('{0}/{1}.fits     {2}\n'.format(path.raw, f, 'IRD_DARK_RAW'))
                    file.close()

                    # products
                    if ctype == 'SKY':
                        loc = 'sky'
                    else:
                        loc = 'internal'
                    dark_file = 'dark_{0}_filt={1}_DIT={2:.2f}'.format(loc, cfilt, DIT)
                    bpm_file  = 'dark_{0}_bpm_filt={1}_DIT={2:.2f}'.format(loc, cfilt, DIT)

                    # different max level in LRS
                    max_level = 1000
                    if cfilt in ['S_LR']:
                        max_level = 15000

                    # esorex parameters
                    args = ['esorex',
                            '--no-checksum=TRUE',
                            '--no-datamd5=TRUE',
                            'sph_ird_master_dark',
                            '--ird.master_dark.sigma_clip=5.0',
                            '--ird.master_dark.save_addprod=TRUE',
                            '--ird.master_dark.max_acceptable={0}'.format(max_level),
                            '--ird.master_dark.outfilename={0}/{1}.fits'.format(path.calib, dark_file),
                            '--ird.master_dark.badpixfilename={0}/{1}.fits'.format(path.calib, bpm_file),
                            sof]

                    # check esorex
                    if shutil.which('esorex') is None:
                        self._logger.error('esorex does not appear to be in your PATH. Please make sure that the ESO pipeline is properly installed before running vlt-sphere.')
                        self._update_recipe_status('sph_ird_cal_dark', sphere.ERROR)
                        return

                    # execute esorex
                    self._logger.debug('> execute esorex')
                    if silent:
                        proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
                    else:
                        proc = subprocess.run(args, cwd=path.tmp)

                    if proc.returncode != 0:
                        self._logger.error('esorex process was not successful')
                        self._update_recipe_status('sph_ird_cal_dark', sphere.ERROR)
                        return

                    # store products
                    self._logger.debug('> update files_info data frame')
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
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')

        # update recipe execution
        self._update_recipe_status('sph_ird_cal_dark', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_cal_detector_flat(self, silent=True):
        '''
        Create the detector flat calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        self._logger.info('Instrument flats')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_cal_detector_flat',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        files_info = self._files_info

        # get list of files
        calibs = files_info[np.logical_not(files_info['PROCESSED']) &
                            (files_info['DPR TYPE'] == 'FLAT,LAMP')]
        filter_combs = calibs['INS COMB IFLT'].unique()

        for cfilt in filter_combs:
            for i, dpr_tech in enumerate(calibs['DPR TECH']):
                if dpr_tech != 'SPECTRUM':
                    date_obs = calibs['DATE-OBS'][i]

                    self._logger.warning(f'The \'DPR TECH\' value of the flat calibration '
                                         f'obtained on {date_obs} is {dpr_tech}. It is '
                                         f'recommended to use flat calibrations for which '
                                         f'the \'DPR TECH\' is \'SPECTRUM\'.')

            cfiles = calibs[calibs['INS COMB IFLT'] == cfilt]
            files = cfiles.index

            self._logger.info(' * filter {0} ({1} files)'.format(cfilt, len(cfiles)))

            # create sof
            self._logger.debug('> create sof file')
            sof = path.sof / 'flat_filt={0}.sof'.format(cfilt)
            file = open(sof, 'w')
            for f in files:
                file.write('{0}/{1}.fits     {2}\n'.format(path.raw, f, 'IRD_FLAT_FIELD_RAW'))
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
                    '--ird.instrument_flat.outfilename={0}/{1}.fits'.format(path.calib, flat_file),
                    '--ird.instrument_flat.badpixfilename={0}/{1}.fits'.format(path.calib, bpm_file),
                    sof]

            # check esorex
            if shutil.which('esorex') is None:
                self._logger.error('esorex does not appear to be in your PATH. Please make sure that the ESO pipeline is properly installed before running vlt-sphere.')
                self._update_recipe_status('sph_ird_cal_detector_flat', sphere.ERROR)
                return

            # execute esorex
            self._logger.debug('> execute esorex')
            if silent:
                proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
            else:
                proc = subprocess.run(args, cwd=path.tmp)

            if proc.returncode != 0:
                self._logger.error('esorex process was not successful')
                self._update_recipe_status('sph_ird_cal_detector_flat', sphere.ERROR)
                return

            # store products
            self._logger.debug('> update files_info data frame')
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
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')

        # update recipe execution
        self._update_recipe_status('sph_ird_cal_detector_flat', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_cal_wave(self, silent=True):
        '''
        Create the wavelength calibration

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        self._logger.info('Wavelength calibration')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_cal_wave',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        files_info = self._files_info

        # get list of files
        wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'LAMP,WAVE')]
        if len(wave_file) != 1:
            self._logger.error('There should be exactly 1 raw wavelength calibration file. Found {0}.'.format(len(wave_file)))
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        DIT = wave_file['DET SEQ1 DIT'][0]
        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_MASTER_DARK') &
                               (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
        if len(dark_file) == 0:
            self._logger.error('There should at least 1 dark file for wavelength calibration. Found none.')
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        filter_comb = wave_file['INS COMB IFLT'][0]
        flat_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_FLAT_FIELD')]
        if len(flat_file) == 0:
            self._logger.error('There should at least 1 flat file for wavelength calibration. Found none.')
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        bpm_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_NON_LINEAR_BADPIXELMAP')]
        if len(flat_file) == 0:
            self._logger.error('There should at least 1 bad pixel map file for wavelength calibration. Found none.')
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        # products
        wav_file = 'wave_calib'

        # laser wavelengths
        wave_lasers = self._wave_cal_lasers

        # esorex parameters
        self._logger.debug('> filter combination is {}'.format(filter_comb))
        if filter_comb == 'S_LR':
            # create standard sof in LRS
            self._logger.debug('> create sof file')
            sof = path.sof / 'wave.sof'
            file = open(sof, 'w')
            file.write('{0}/{1}.fits     {2}\n'.format(path.raw, wave_file.index[0], 'IRD_WAVECALIB_RAW'))
            file.write('{0}/{1}.fits     {2}\n'.format(path.calib, dark_file.index[0], 'IRD_MASTER_DARK'))
            file.write('{0}/{1}.fits     {2}\n'.format(path.calib, flat_file.index[0], 'IRD_FLAT_FIELD'))
            file.write('{0}/{1}.fits     {2}\n'.format(path.calib, bpm_file.index[0], 'IRD_STATIC_BADPIXELMAP'))
            file.close()

            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_wave_calib',
                    '--ird.wave_calib.column_width=200',
                    '--ird.wave_calib.grism_mode=FALSE',
                    '--ird.wave_calib.threshold=2000',
                    '--ird.wave_calib.number_lines=6',
                    '--ird.wave_calib.wavelength_line1={:.2f}'.format(wave_lasers[0]),
                    '--ird.wave_calib.wavelength_line2={:.2f}'.format(wave_lasers[1]),
                    '--ird.wave_calib.wavelength_line3={:.2f}'.format(wave_lasers[2]),
                    '--ird.wave_calib.wavelength_line4={:.2f}'.format(wave_lasers[3]),
                    '--ird.wave_calib.wavelength_line5={:.2f}'.format(wave_lasers[4]),
                    '--ird.wave_calib.wavelength_line6={:.2f}'.format(wave_lasers[5]),
                    '--ird.wave_calib.outfilename={0}/{1}.fits'.format(path.calib, wav_file),
                    sof]
        elif filter_comb == 'S_MR':
            # masking of second order spectrum in MRS
            self._logger.debug('> masking second order')
            wave_fname = wave_file.index[0]
            wave_data, hdr = fits.getdata(path.raw / '{}.fits'.format(wave_fname), header=True)
            wave_data = wave_data.squeeze()
            wave_data[:60, :] = 0
            fits.writeto(path.preproc / '{}_masked.fits'.format(wave_fname), wave_data, hdr, overwrite=True,
                         output_verify='silentfix')

            # create sof using the masked file
            self._logger.debug('> create sof file')
            sof = path.sof / 'wave.sof'
            file = open(sof, 'w')
            file.write('{0}/{1}_masked.fits {2}\n'.format(path.preproc, wave_fname, 'IRD_WAVECALIB_RAW'))
            file.write('{0}/{1}.fits        {2}\n'.format(path.calib, dark_file.index[0], 'IRD_MASTER_DARK'))
            file.write('{0}/{1}.fits        {2}\n'.format(path.calib, flat_file.index[0], 'IRD_FLAT_FIELD'))
            file.write('{0}/{1}.fits        {2}\n'.format(path.calib, bpm_file.index[0], 'IRD_STATIC_BADPIXELMAP'))
            file.close()

            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_wave_calib',
                    '--ird.wave_calib.column_width=200',
                    '--ird.wave_calib.grism_mode=TRUE',
                    '--ird.wave_calib.threshold=1000',
                    '--ird.wave_calib.number_lines=5',
                    '--ird.wave_calib.wavelength_line1={:.2f}'.format(wave_lasers[0]),
                    '--ird.wave_calib.wavelength_line2={:.2f}'.format(wave_lasers[1]),
                    '--ird.wave_calib.wavelength_line3={:.2f}'.format(wave_lasers[2]),
                    '--ird.wave_calib.wavelength_line4={:.2f}'.format(wave_lasers[3]),
                    '--ird.wave_calib.wavelength_line5={:.2f}'.format(wave_lasers[4]),
                    '--ird.wave_calib.outfilename={0}/{1}.fits'.format(path.calib, wav_file),
                    sof]

        # check esorex
        if shutil.which('esorex') is None:
            self._logger.error('esorex does not appear to be in your PATH. Please make sure that the ESO pipeline is properly installed before running vlt-sphere.')
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        # execute esorex
        self._logger.debug('> execute esorex')
        if silent:
            proc = subprocess.run(args, cwd=path.tmp, stdout=subprocess.DEVNULL)
        else:
            proc = subprocess.run(args, cwd=path.tmp)

        if proc.returncode != 0:
            self._logger.error('esorex process was not successful')
            self._update_recipe_status('sph_ird_cal_wave', sphere.ERROR)
            return

        # store products
        self._logger.debug('> update files_info data frame')
        files_info.loc[wav_file, 'DPR CATG'] = wave_file['DPR CATG'][0]
        files_info.loc[wav_file, 'DPR TYPE'] = wave_file['DPR TYPE'][0]
        files_info.loc[wav_file, 'INS COMB IFLT'] = wave_file['INS COMB IFLT'][0]
        files_info.loc[wav_file, 'INS4 FILT2 NAME'] = wave_file['INS4 FILT2 NAME'][0]
        files_info.loc[wav_file, 'INS1 MODE'] = wave_file['INS1 MODE'][0]
        files_info.loc[wav_file, 'INS1 FILT NAME'] = wave_file['INS1 FILT NAME'][0]
        files_info.loc[wav_file, 'INS1 OPTI2 NAME'] = wave_file['INS1 OPTI2 NAME'][0]
        files_info.loc[wav_file, 'DET SEQ1 DIT'] = wave_file['DET SEQ1 DIT'][0]
        files_info.loc[wav_file, 'PROCESSED'] = True
        files_info.loc[wav_file, 'PRO CATG'] = 'IRD_WAVECALIB'

        # save
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')

        # store default wavelength calibration in preproc
        self._logger.debug('> compute default wavelength calibration')
        if filter_comb == 'S_LR':
            centers  = self._default_center_lrs
            wave_min = self._wave_min_lrs
            wave_max = self._wave_max_lrs
        elif filter_comb == 'S_MR':
            centers  = self._default_center_mrs
            wave_min = self._wave_min_mrs
            wave_max = self._wave_max_mrs

        wave_calib = fits.getdata(path.calib / '{}.fits'.format(wav_file))
        wave_lin   = get_wavelength_calibration(wave_calib, centers, wave_min, wave_max)

        self._logger.debug('> save default wavelength calibration')
        fits.writeto(path.preproc / 'wavelength_default.fits', wave_lin.T, overwrite=True)

        # update recipe execution
        self._update_recipe_status('sph_ird_cal_wave', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_preprocess_science(self,
                                   subtract_background=True, fix_badpix=True,
                                   collapse_science=False, collapse_psf=True, collapse_center=True):
        '''Pre-processes the science frames.

        This function can perform multiple steps:
          - collapse of the frames
          - subtract the background
          - correct bad pixels
          - reformat IRDIS data in (x,y,lambda) cubes

        For the science, PSFs or star center frames, the full cubes
        are mean-combined into a single frame.

        The pre-processed frames are saved in the preproc
        sub-directory and will be combined later.

        Parameters
        ----------
        subtract_background : bool
            Performs background subtraction. Default is True

        fix_badpix : bool
            Performs correction of bad pixels. Default is True

        collapse_science :  bool
            Collapse data for OBJECT cubes. Default is False

        collapse_psf :  bool
            Collapse data for OBJECT,FLUX cubes. Default is True. Note
            that the collapse type is mean and cannot be changed.

        collapse_center :  bool
            Collapse data for OBJECT,CENTER cubes. Default is True. Note
            that the collapse type is mean and cannot be changed.

        '''

        self._logger.info('Pre-process science files')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_preprocess_science',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        files_info = self._files_info
        frames_info = self._frames_info

        # clean before we start
        self._logger.debug('> remove old preproc files')
        files = path.preproc.glob('*_DIT???_preproc.fits')
        for file in files:
            file.unlink()

        # filter combination
        filter_comb = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()[0]

        # bpm
        if fix_badpix:
            bpm_files = files_info[(files_info['PRO CATG'] == 'IRD_STATIC_BADPIXELMAP') |
                                   (files_info['PRO CATG'] == 'IRD_NON_LINEAR_BADPIXELMAP')].index
            bpm_files = [path.calib / '{}.fits'.format(f) for f in bpm_files]
            if len(bpm_files) == 0:
                self._logger.error('Could not fin any bad pixel maps')
                self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                return

            bpm = toolbox.compute_bad_pixel_map(bpm_files, logger=self._logger)

            # mask dead regions
            bpm[:15, :]      = 0
            bpm[1013:, :]    = 0
            bpm[:, :50]      = 0
            bpm[:, 941:1078] = 0
            bpm[:, 1966:]    = 0

        # flat
        flat_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_FLAT_FIELD') &
                               (files_info['INS COMB IFLT'] == filter_comb)]
        if len(flat_file) != 1:
            self._logger.error('There should be exactly 1 flat file. Found {0}.'.format(len(flat_file)))
            self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
            return
        flat = fits.getdata(path.calib / '{}.fits'.format(flat_file.index[0]))

        # final dataframe
        self._logger.debug('> create frames_info_preproc data frame')
        index = pd.MultiIndex(names=['FILE', 'IMG'], levels=[[], []], codes=[[], []])
        frames_info_preproc = pd.DataFrame(index=index, columns=frames_info.columns, dtype='float')

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

                self._logger.info('{0} files of type {1} with DIT={2} sec'.format(len(sfiles), typ, DIT))

                if subtract_background:
                    # look for sky, then background, then darks
                    # normally there should be only one with the proper DIT
                    dfiles = []
                    for d in dark_types:
                        dfiles = files_info[(files_info['PRO CATG'] == 'IRD_MASTER_DARK') &
                                            (files_info['DPR TYPE'] == d) & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
                        if len(dfiles) != 0:
                            break
                    self._logger.info('   ==> found {0} corresponding {1} file'.format(len(dfiles), d))

                    if len(dfiles) == 0:
                        # issue a warning if absolutely no background is found
                        self._logger.warning('No background has been found. Pre-processing will continue but data quality will likely be affected')
                        bkg = np.zeros((1024, 2048))
                    elif len(dfiles) == 1:
                        bkg = fits.getdata(path.calib / '{}.fits'.format(dfiles.index[0]))
                    elif len(dfiles) > 1:
                        # FIXME: handle cases when multiple backgrounds are found?
                        self._logger.error('Unexpected number of background files ({0})'.format(len(dfiles)))
                        self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                        return

                # process files
                for idx, (fname, finfo) in enumerate(sfiles.iterrows()):
                    # frames_info extract
                    finfo = frames_info.loc[(fname, slice(None)), :]

                    self._logger.info(' * file {0}/{1}: {2}, NDIT={3}'.format(idx+1, len(sfiles), fname, len(finfo)))

                    # read data
                    self._logger.info('   ==> read data')
                    img, hdr = fits.getdata(path.raw / '{}.fits'.format(fname), header=True)

                    # add extra dimension to single images to make cubes
                    if img.ndim == 2:
                        img = img[np.newaxis, ...]

                    # mask dead regions
                    img[:, :15, :]      = np.nan
                    img[:, 1013:, :]    = np.nan
                    img[:, :, :50]      = np.nan
                    img[:, :, 941:1078] = np.nan
                    img[:, :, 1966:]    = np.nan

                    # collapse
                    if (typ == 'OBJECT,CENTER'):
                        if collapse_center:
                            self._logger.info('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean', logger=self._logger)
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none', logger=self._logger)
                    elif (typ == 'OBJECT,FLUX'):
                        if collapse_psf:
                            self._logger.info('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean', logger=self._logger)
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none', logger=self._logger)
                    elif (typ == 'OBJECT'):
                        if collapse_science:
                            self._logger.info('   ==> collapse: mean ({0} -> 1 frame, 0 dropped)'.format(len(img)))
                            img = np.mean(img, axis=0, keepdims=True)

                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'mean', logger=self._logger)
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none', logger=self._logger)

                    # check for any error during collapse of frame information
                    if frames_info_new is None:
                        self._logger.error('An error occured when collapsing frames info')
                        self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                        return
                    
                    # merge frames info
                    frames_info_preproc = pd.concat((frames_info_preproc, frames_info_new))

                    # background subtraction
                    if subtract_background:
                        self._logger.info('   ==> subtract background')
                        for f in range(len(img)):
                            img[f] -= bkg

                    # divide flat
                    if subtract_background:
                        self._logger.info('   ==> divide by flat field')
                        for f in range(len(img)):
                            img[f] /= flat

                    # bad pixels correction
                    if fix_badpix:
                        self._logger.info('   ==> correct bad pixels')
                        for f in range(len(img)):
                            frame = img[f]
                            frame = imutils.fix_badpix(frame, bpm, npix=12, weight=True)

                            # additional sigma clipping to remove transitory bad pixels
                            # not done for OBJECT,FLUX because PSF peak can be clipped
                            if (typ != 'OBJECT,FLUX'):
                                frame = imutils.sigma_filter(frame, box=7, nsigma=4, iterate=False)

                            img[f] = frame

                    # reshape data
                    self._logger.info('   ==> reshape data')
                    NDIT = img.shape[0]
                    nimg = np.zeros((NDIT, 2, 1024, 1024))
                    for f in range(len(img)):
                        nimg[f, 0] = img[f, :, 0:1024]
                        nimg[f, 1] = img[f, :, 1024:]
                    img = nimg

                    # save DITs individually
                    self._logger.debug('> save pre-processed images')
                    for f in range(len(img)):
                        frame = nimg[f, ...].squeeze()
                        hdr['HIERARCH ESO DET NDIT'] = 1
                        fits.writeto(path.preproc / '{}_DIT{:03d}_preproc.fits'.format(fname, f), frame, hdr,
                                     overwrite=True, output_verify='silentfix')

        # sort and save final dataframe
        self._logger.debug('> save frames_info_preproc.csv')
        frames_info_preproc.sort_values(by='TIME', inplace=True)
        frames_info_preproc.to_csv(path.preproc / 'frames_preproc.csv')

        self._frames_info_preproc = frames_info_preproc

        # update recipe execution
        self._update_recipe_status('sph_ird_preprocess_science', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_star_center(self, high_pass=False, plot=True):
        '''Determines the star center for all frames where a center can be
        determined (OBJECT,CENTER and OBJECT,FLUX)

        Parameters
        ----------
        high_pass : bool
            Apply high-pass filter to the image before searching for the satelitte spots.
            Default is False

        plot : bool
            Display and save diagnostic plot for quality check. Default is True

        '''

        self._logger.info('Star centers determination')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_star_center',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        pixel = self._pixel
        files_info  = self._files_info
        frames_info = self._frames_info_preproc

        # resolution-specific parameters
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        if filter_comb == 'S_LR':
            centers  = self._default_center_lrs
            wave_min = self._wave_min_lrs
            wave_max = self._wave_max_lrs
        elif filter_comb == 'S_MR':
            centers  = self._default_center_mrs
            wave_min = self._wave_min_mrs
            wave_max = self._wave_max_mrs

        # wavelength map
        self._logger.debug('> compute default wavelength calibration')
        wave_file  = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_WAVECALIB')]
        wave_calib = fits.getdata(path.calib / '{}.fits'.format(wave_file.index[0]))
        wave_lin   = get_wavelength_calibration(wave_calib, centers, wave_min, wave_max)

        # start with OBJECT,FLUX
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        if len(flux_files) != 0:
            for file, idx in flux_files.index:
                self._logger.info(' * OBJECT,FLUX: {0}'.format(file))

                # read data
                self._logger.debug('> read data')
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube, hdr = fits.getdata(path.preproc / '{}.fits'.format(fname), header=True)

                # centers
                if plot:
                    save_path = path.products / '{}_PSF_fitting.pdf'.format(fname)
                else:
                    save_path = None
                psf_center = toolbox.star_centers_from_PSF_lss_cube(cube, wave_lin, pixel, save_path=save_path,
                                                                    logger=self._logger)

                # save
                self._logger.debug('> save centers')
                fits.writeto(path.preproc / '{}_centers.fits'.format(fname), psf_center, overwrite=True)

        # then OBJECT,CENTER (if any)
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) != 0:
            DIT = starcen_files['DET SEQ1 DIT'].round(2)[0]
            starsci_files = frames_info[(frames_info['DPR TYPE'] == 'OBJECT') & (frames_info['DET SEQ1 DIT'].round(2) == DIT)]

            for file, idx in starcen_files.index:
                self._logger.info(' * OBJECT,CENTER: {0}'.format(file))

                # read center data
                self._logger.debug('> read data')
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube_cen, hdr = fits.getdata(path.preproc / '{}.fits'.format(fname), header=True)

                # read science data
                self._logger.debug('> read matching science data')
                if len(starsci_files) != 0:
                    self._logger.debug('> read matiching science data')
                    fname2 = '{0}_DIT{1:03d}_preproc'.format(starsci_files.index[0][0], idx)
                    cube_sci, hdr = fits.getdata(path.preproc / '{}.fits'.format(fname2), header=True)
                else:
                    cube_sci = None

                # centers
                if plot:
                    save_path = path.products / '{}_spots_fitting.pdf'.format(fname)
                else:
                    save_path = None
                spot_centers, spot_dist, img_centers \
                    = toolbox.star_centers_from_waffle_lss_cube(cube_cen, cube_sci, wave_lin, centers, pixel,
                                                                high_pass=high_pass, save_path=save_path,
                                                                logger=self._logger)

                # save
                self._logger.debug('> save centers')
                fits.writeto(path.preproc / '{}_centers.fits'.format(fname), img_centers, overwrite=True)
                fits.writeto(path.preproc / '{}_spot_distance.fits'.format(fname), spot_dist, overwrite=True)

        # update recipe execution
        self._update_recipe_status('sph_ird_star_center', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE

        
    def sph_ird_wavelength_recalibration(self, fit_scaling=True, plot=True):
        '''Performs a recalibration of the wavelength, if star center frames
        are available. Otherwise simply use the ESO pipeline-calibrated law.

        It follows a similar process to that used for the IFS
        data. The method for the IFS is described in Vigan et
        al. (2015, MNRAS, 454, 129):

        https://ui.adsabs.harvard.edu/#abs/2015MNRAS.454..129V/abstract

        Parameters
        ----------
        fit_scaling : bool
            Perform a polynomial fitting of the wavelength scaling
            law. It helps removing high-frequency noise that can
            result from the waffle fitting. Default is True

        plot : bool
            Display and save diagnostic plot for quality check. Default is True

        '''

        self._logger.info('Wavelength recalibration')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_wavelength_recalibration',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        lasers = self._wave_cal_lasers
        files_info  = self._files_info
        frames_info = self._frames_info_preproc

        # remove old files
        self._logger.debug('> remove old recalibrated wavelength calibration')
        wfile = path.preproc / 'wavelength_recalibrated.fits'
        if wfile.exists():
            wfile.unlink()

        # resolution-specific parameters
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        if filter_comb == 'S_LR':
            centers  = self._default_center_lrs
            wave_min = self._wave_min_lrs
            wave_max = self._wave_max_lrs
        elif filter_comb == 'S_MR':
            centers  = self._default_center_mrs
            wave_min = self._wave_min_mrs
            wave_max = self._wave_max_mrs

        # wavelength map
        self._logger.debug('> compute default wavelength calibration')
        wave_file  = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_WAVECALIB')]
        wave_calib = fits.getdata(path.calib / '{}.fits'.format(wave_file.index[0]))
        wave_lin   = get_wavelength_calibration(wave_calib, centers, wave_min, wave_max)

        # reference wavelength
        idx_ref = 3
        wave_ref = lasers[idx_ref]

        # get spot distance from the first OBJECT,CENTER in the sequence
        self._logger.debug('> read saved spot distances')
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) == 0:
            self._logger.info('   ==> no OBJECT,CENTER file in the data set. Wavelength cannot be recalibrated. The standard wavelength calibrated by the ESO pripeline will be used.')
            return

        fname = '{0}_DIT{1:03d}_preproc_spot_distance'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
        spot_dist = fits.getdata(path.preproc / '{}.fits'.format(fname))

        if plot:
            pdf = PdfPages(path.products / 'wavelength_recalibration.pdf')

        pix = np.arange(1024)
        wave_final = np.zeros((1024, 2))
        for fidx in range(2):
            self._logger.info('   ==> field {0:2d}/{1:2d}'.format(fidx+1, 2))

            wave = wave_lin[fidx]
            dist = spot_dist[:, fidx]

            imin = np.nanargmin(np.abs(wave-wave_ref))

            # scaling factor
            self._logger.debug('> compute wavelength scaling factor')
            scaling_raw = dist / dist[imin]

            self._logger.debug('> polynomial fit for recalibration')
            if filter_comb == 'S_LR':
                # FIXME: implement smoothing of the scaling factor for
                # LRS mode
                self._logger.error('Wavelength recalibration is not yet implemented for IRDIS-LRS mode')
                self._update_recipe_status('sph_ird_wavelength_recalibration', sphere.ERROR)
                return
            elif filter_comb == 'S_MR':
                # linear fit with a 5-degree polynomial
                good = np.where(np.isfinite(wave))
                p = np.polyfit(pix[good], scaling_raw[good], 5)

                scaling_fit = np.polyval(p, pix)

            wave_final_raw = wave[imin] * scaling_raw
            wave_final_fit = wave[imin] * scaling_fit

            bad = np.where(np.logical_not(np.isfinite(wave)))
            wave_final_raw[bad] = np.nan
            wave_final_fit[bad] = np.nan

            wave_diff = np.abs(wave_final_fit - wave)
            self._logger.info('   ==> difference with calibrated wavelength: min={0:.1f} nm, max={1:.1f} nm'.format(np.nanmin(wave_diff), np.nanmax(wave_diff)))

            if fit_scaling:
                self._logger.debug('> use fitted scaling factor')
                wave_final[:, fidx] = wave_final_fit
                use_r = ''
                use_f = ' <=='
            else:
                self._logger.debug('> use raw scaling factor')
                wave_final[:, fidx] = wave_final_raw
                use_r = ' <=='
                use_f = ''

            # plot
            if plot:
                plt.figure('Wavelength recalibration', figsize=(10, 10))
                plt.clf()

                plt.subplot(211)
                plt.axvline(imin, color='k', linestyle='--')
                plt.plot(pix, wave, label='DRH', color='r', lw=3)
                plt.plot(pix, wave_final_raw, label='Recalibrated [raw]'+use_r)
                plt.plot(pix, wave_final_fit, label='Recalibrated [fit]'+use_f)
                plt.legend(loc='upper left')
                plt.ylabel('Wavelength r[nm]')
                plt.title('Field #{}'.format(fidx))
                plt.xlim(1024, 0)
                plt.gca().xaxis.set_ticklabels([])

                plt.subplot(212)
                plt.axvline(imin, color='k', linestyle='--')
                plt.plot(pix, wave-wave_final_raw)
                plt.plot(pix, wave-wave_final_fit)
                plt.ylabel('Residuals r[nm]')
                plt.xlabel('Detector coordinate [pix]')
                plt.xlim(1024, 0)

                plt.subplots_adjust(left=0.13, right=0.97, bottom=0.08, top=0.96, hspace=0.05)

                pdf.savefig()

        if plot:
            pdf.close()

        # save
        self._logger.info(' * saving')
        fits.writeto(path.preproc / 'wavelength_recalibrated.fits', wave_final, overwrite=True)

        # update recipe execution
        self._update_recipe_status('sph_ird_wavelength_recalibration', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_combine_data(self, cpix=True, psf_dim=80, science_dim=800, correct_mrs_chromatism=True,
                             split_posang=True, shift_method='fft', manual_center=None, coarse_centering=False):
        '''Combine and save the science data into final cubes

        All types of data are combined independently: PSFs
        (OBJECT,FLUX), star centers (OBJECT,CENTER) and standard
        coronagraphic images (OBJECT).

        Depending on the observing strategy, there can be several
        position angle positions in the sequence. Images taken at
        different position angles can be either kept together or
        split into different cubes. In either case a posang vector
        is saved alongside the science cube(s).

        For each type of data, the method saves 3 different files:

          - *_cube: the (x,y,time,nfield) cube

          - *_posang: the position angle vector.

          - *_frames: a csv file with all the information for every
                      frames. There is one line by time step in the
                      data cube.

        Centering
        ---------

        By default, a fine (sub-pixel) centering is performed if the
        an OBJECT,CENTER frame was acquired in the sequence or if
        there is a valid user-provided center. However, if the
        coarse_centering keyword is set to True, only a "coarse
        centering" is performed, which requires no interpolation:

          - only integer shifts (shift_method='roll')
          - centering on an integer pixel (cpix=True)
          - no correction of the MRS chromatism (correct_mrs_chromatism=False)

        This option is useful if the user wants to perform a
        posteriori centering of the spectrum, e.g. to fully preserve
        photometry.

        If there was no OBJECT,CENTER acquired in the sequence, then
        the centering will be performed with respect to a default,
        pre-defined center that is representative of the typical center
        of the coronagraph.

        Parameters
        ----------
        cpix : bool
            If True the images are centered on the pixel at coordinate
            dim//2 in the spatial dimension. If False the images are
            centered between 2 pixels, at coordinates (dim-1)/2. The
            value of cpix is automatically set to True when
            coarse_centering is set to True. Default is True.

        psf_dim : even int
            Size of the PSF images along in the spatial
            dimension. Default is 80x pixels

        science_dim : even int
            Size of the science images (star centers and standard
            coronagraphic images) in the spatial dimension. Default is
            800 pixels

        correct_mrs_chromatism : bool
            Correct for the slight chromatism in the MRS mode. This
            chromatism induces a slight shift of the PSF center with
            wavelength. The value of correct_mrs_chromatism is
            automatically set to True when coarse_centering is set to
            True. Default is True.

        split_posang : bool
            Save data taken at different position angles in separate
            science files. Default is True

        manual_center : array
            User provided spatial center for the OBJECT,CENTER and
            OBJECT frames. This should be an array of 2 values (cx for
            the 2 IRDIS fields). Default is None

        coarse_centering : bool
            Control if images are finely centered or not before being
            combined. However the images are still roughly centered by
            shifting them by an integer number of pixel to bring the
            center of the data close to the center of the images. This
            option is useful if fine centering must be done
            afterwards. Default is False.

        shift_method : str
            Method to shifting the images: fft or interp.  Default is
            fft

        '''

        self._logger.info('Combine science data')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_combine_data',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path
        nwave = self._nwave
        frames_info = self._frames_info_preproc

        # resolution-specific parameters
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        if filter_comb == 'S_LR':
            default_center = self._default_center_lrs
            wave_min = self._wave_min_lrs
            wave_max = self._wave_max_lrs
        elif filter_comb == 'S_MR':
            default_center = self._default_center_mrs
            wave_min = self._wave_min_mrs
            wave_max = self._wave_max_mrs

        # read final wavelength calibration
        self._logger.debug('> save final wavelength')
        wfile = path.preproc / 'wavelength_recalibrated.fits'
        if wfile.exists():
            wave = fits.getdata(wfile)
        else:
            wfile = path.preproc / 'wavelength_default.fits'
            if wfile.exists():
                self._logger.warning('Using default wavelength calibration.')
                wave = fits.getdata(wfile)
            else:
                self._logger.error('Missing default or recalibrated wavelength calibration. You must first run either sph_ird_cal_wave or sph_ird_wavelength_recalibration().')
                self._update_recipe_status('sph_ird_combine_data', sphere.ERROR)
                return

        # wavelength solution: make sure we have the same number of
        # wave points in each field
        mask   = ((wave_min <= wave) & (wave <= wave_max))
        iwave0 = np.where(mask[:, 0])[0]
        iwave1 = np.where(mask[:, 1])[0]
        nwave  = np.min([iwave0.size, iwave1.size])

        iwave = np.empty((nwave, 2), dtype=np.int)
        iwave[:, 0] = iwave0[:nwave]
        iwave[:, 1] = iwave1[:nwave]

        final_wave = np.empty((nwave, 2))
        final_wave[:, 0] = wave[iwave[:, 0], 0]
        final_wave[:, 1] = wave[iwave[:, 1], 1]

        fits.writeto(path.products / 'wavelength.fits', final_wave.squeeze().T, overwrite=True)

        # max images size
        if psf_dim > 1024:
            self._logger.warning('psf_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            psf_dim = 1024

        if science_dim > 1024:
            self._logger.warning('science_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            science_dim = 1024

        # centering configuration
        if coarse_centering:
            self._logger.warning('Images will be coarsely centered without any interpolation. Automatic settings for coarse centering: shift_method=\'roll\', cpix=True, correct_mrs_chromatism=False')
            shift_method = 'roll'
            cpix = True
            correct_mrs_chromatism = False

        if manual_center is not None:
            manual_center = np.array(manual_center)

            if manual_center.shape != (2,):
                self._logger.error('manual_center does not have the right number of dimensions.')
                self._update_recipe_status('sph_ird_combine_data', sphere.ERROR)
                return

            self._logger.warning('Images will be centered using the user-provided center ({},{})'.format(*manual_center))

            manual_center = np.full((1024, 2), manual_center, dtype=np.float)

        #
        # OBJECT,FLUX
        #
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        nfiles = len(flux_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT,FLUX data')

            # final arrays
            psf_cube   = np.zeros((2, nfiles, nwave, psf_dim))
            psf_posang = np.zeros(nfiles)

            # final center
            if cpix:
                cc = psf_dim // 2
            else:
                cc = (psf_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(flux_files.index):
                self._logger.info('   ==> file {0}/{1}: {2}, DIT #{3}'.format(file_idx+1, len(flux_files), file, idx))

                # read data
                self._logger.debug('> read data')
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube = fits.getdata(path.preproc / '{}.fits'.format(fname))

                self._logger.debug('> read centers')
                cfile = path.preproc / '{}_centers.fits'.format(fname)
                if cfile.exists():
                    centers = fits.getdata(cfile)
                else:
                    self._logger.warning('sph_ird_star_center() has not been executed. Images will be centered using default centers ({}, {})'.format(*default_center[:, 0]))
                    centers = np.full((1024, 2), default_center[:, 0], dtype=np.float)

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                psf_posang[file_idx] = frames_info.loc[(file, idx), 'INS4 DROT2 POSANG'] + 90

                # center
                for field_idx, img in enumerate(cube):
                    self._logger.debug('> field {}'.format(field_idx))
                    # wavelength solution for this field
                    ciwave = iwave[:, field_idx]

                    if correct_mrs_chromatism and (filter_comb == 'S_MR'):
                        self._logger.debug('> correct MRS chromatism')
                        img = img.astype(np.float)
                        self._logger.debug('> shift and normalize')
                        for wave_idx, widx in enumerate(ciwave):
                            cx = centers[widx, field_idx]

                            line = img[widx, :]

                            nimg = imutils.shift(line, cc-cx, method=shift_method)
                            nimg = nimg / DIT

                            psf_cube[field_idx, file_idx, wave_idx] = nimg[:psf_dim]
                    else:
                        cx = centers[ciwave, field_idx].mean()

                        self._logger.debug('> shift and normalize')
                        img  = img.astype(np.float)
                        nimg = imutils.shift(img, (cc-cx, 0), method=shift_method)
                        nimg = nimg / DIT

                        psf_cube[field_idx, file_idx] = nimg[ciwave, :psf_dim]

                    # neutral density
                    self._logger.debug('> compensate for neutral density')
                    cwave  = final_wave[:, field_idx]
                    ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                    w, attenuation = transmission.transmission_nd(ND, wave=cwave)
                    psf_cube[field_idx, file_idx] = (psf_cube[field_idx, file_idx].T / attenuation).T

            self._logger.debug('> save final cubes and metadata')
            if split_posang:
                self._logger.debug('> split position angles')
                pas = np.unique(psf_posang)
                for pa in pas:
                    ii = np.where(psf_posang == pa)[0]

                    # save metadata
                    flux_files[(flux_files['INS4 DROT2 POSANG'] + 90) == pa].to_csv(path.products / 'psf_posang={:06.2f}_frames.csv'.format(pa))
                    fits.writeto(path.products / 'psf_posang={:06.2f}_posang.fits'.format(pa), psf_posang[ii], overwrite=True)

                    # save final cubes
                    fits.writeto(path.products / 'psf_posang={:06.2f}_cube.fits'.format(pa), psf_cube[:, ii], overwrite=True)
            else:
                # save metadata
                flux_files.to_csv(path.products / 'psf_posang=all_frames.csv')
                fits.writeto(path.products / 'psf_posang=all_posang.fits', psf_posang, overwrite=True)

                # save final cubes
                fits.writeto(path.products / 'psf_posang=all_cube.fits', psf_cube, overwrite=True)

            # delete big cubes
            self._logger.debug('> free memory')
            del psf_cube

        #
        # OBJECT,CENTER
        #
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        nfiles = len(starcen_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT,CENTER data')

            # final arrays
            cen_cube   = np.zeros((2, nfiles, nwave, science_dim))
            cen_posang = np.zeros(nfiles)

            # final center
            if cpix:
                cc = science_dim // 2
            else:
                cc = (science_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(starcen_files.index):
                self._logger.info('   ==> file {0}/{1}: {2}, DIT #{3}'.format(file_idx+1, len(starcen_files), file, idx))

                # read data
                self._logger.debug('> read data')
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube = fits.getdata(path.preproc / '{}.fits'.format(fname))

                # use manual center if explicitely requested
                self._logger.debug('> read centers')
                if manual_center is not None:
                    centers = manual_center
                else:
                    centers = fits.getdata(path.preproc / '{}_centers.fits'.format(fname))

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                cen_posang[file_idx] = frames_info.loc[(file, idx), 'INS4 DROT2 POSANG'] + 90

                # center
                for field_idx, img in enumerate(cube):
                    self._logger.debug('> field {}'.format(field_idx))
                    # wavelength solution for this field
                    ciwave = iwave[:, field_idx]

                    if correct_mrs_chromatism and (filter_comb == 'S_MR'):
                        self._logger.debug('> correct MRS chromatism')
                        img = img.astype(np.float)
                        self._logger.debug('> shift and normalize')
                        for wave_idx, widx in enumerate(ciwave):
                            cx = centers[widx, field_idx]

                            line = img[widx, :]
                            nimg = imutils.shift(line, cc-cx, method=shift_method)
                            nimg = nimg / DIT

                            cen_cube[field_idx, file_idx, wave_idx] = nimg[:science_dim]
                    else:
                        cx = centers[ciwave, field_idx].mean()

                        self._logger.debug('> shift and normalize')
                        img  = img.astype(np.float)
                        nimg = imutils.shift(img, (cc-cx, 0), method=shift_method)
                        nimg = nimg / DIT

                        cen_cube[field_idx, file_idx] = nimg[ciwave, :science_dim]

                    # neutral density
                    self._logger.debug('> compensate for neutral density')
                    cwave  = final_wave[:, field_idx]
                    ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                    w, attenuation = transmission.transmission_nd(ND, wave=cwave)
                    cen_cube[field_idx, file_idx] = (cen_cube[field_idx, file_idx].T / attenuation).T

            self._logger.debug('> save final cubes and metadata')
            if split_posang:
                self._logger.debug('> split position angles')
                pas = np.unique(cen_posang)
                for pa in pas:
                    ii = np.where(cen_posang == pa)[0]

                    # save metadata
                    starcen_files[(starcen_files['INS4 DROT2 POSANG'] + 90) == pa].to_csv(path.products / 'starcenter_posang={:06.2f}_frames.csv'.format(pa))
                    fits.writeto(path.products / 'starcenter_posang={:06.2f}_posang.fits'.format(pa), cen_posang[ii], overwrite=True)

                    # save final cubes
                    fits.writeto(path.products / 'starcenter_posang={:06.2f}_cube.fits'.format(pa), cen_cube[:, ii], overwrite=True)
            else:
                # save metadata
                starcen_files.to_csv(path.products / 'starcenter_posang=all_frames.csv')
                fits.writeto(path.products / 'starcenter_posang=all_posang.fits', cen_posang, overwrite=True)

                # save final cubes
                fits.writeto(path.products / 'starcenter_posang=all_cube.fits', cen_cube, overwrite=True)

            # delete big cubes
            del cen_cube

        #
        # OBJECT
        #
        object_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT']
        nfiles = len(object_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT data')

            # final arrays
            sci_cube   = np.zeros((2, nfiles, nwave, science_dim))
            sci_posang = np.zeros(nfiles)

            # use manual center if explicitely requested
            self._logger.debug('> read centers')
            if manual_center is not None:
                centers = np.full((1024, 2), manual_center, dtype=np.float)
            else:
                # FIXME: ticket #12. Use first DIT of first OBJECT,CENTER
                # in the sequence, but it would be better to be able to
                # select which CENTER to use
                starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
                if len(starcen_files) == 0:
                    self._logger.warning('No OBJECT,CENTER file in the data set. Images will be centered using default center ({},{})'.format(*default_center[:, 0]))
                    centers = np.full((1024, 2), default_center[:, 0], dtype=np.float)
                else:
                    fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                    centers = fits.getdata(path.preproc / fname)

            # make sure we have only integers if user wants coarse centering
            if coarse_centering:
                centers = centers.astype(np.int)

            # final center
            if cpix:
                cc = science_dim // 2
            else:
                cc = (science_dim - 1) / 2

            # read and combine files
            for file_idx, (file, idx) in enumerate(object_files.index):
                self._logger.info('   ==> file {0}/{1}: {2}, DIT #{3}'.format(file_idx+1, len(object_files), file, idx))

                # read data
                self._logger.debug('> read data')
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube = fits.getdata(path.preproc / '{}.fits'.format(fname))

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                sci_posang[file_idx] = frames_info.loc[(file, idx), 'INS4 DROT2 POSANG'] + 90

                # center
                for field_idx, img in enumerate(cube):
                    self._logger.debug('> field {}'.format(field_idx))
                    # wavelength solution for this field
                    ciwave = iwave[:, field_idx]

                    if correct_mrs_chromatism and (filter_comb == 'S_MR'):
                        self._logger.debug('> correct MRS chromatism')
                        img = img.astype(np.float)
                        self._logger.debug('> shift and normalize')
                        for wave_idx, widx in enumerate(ciwave):
                            cx = centers[widx, field_idx]

                            line = img[widx, :]
                            nimg = imutils.shift(line, cc-cx, method=shift_method)
                            nimg = nimg / DIT

                            sci_cube[field_idx, file_idx, wave_idx] = nimg[:science_dim]
                    else:
                        cx = centers[ciwave, field_idx].mean()

                        self._logger.debug('> shift and normalize')
                        img  = img.astype(np.float)
                        nimg = imutils.shift(img, (cc-cx, 0), method=shift_method)
                        nimg = nimg / DIT

                        sci_cube[field_idx, file_idx] = nimg[ciwave, :science_dim]

                    # neutral density
                    self._logger.debug('> compensate for neutral density')
                    cwave  = final_wave[:, field_idx]
                    ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                    w, attenuation = transmission.transmission_nd(ND, wave=cwave)
                    sci_cube[field_idx, file_idx] = (sci_cube[field_idx, file_idx].T / attenuation).T

            self._logger.debug('> save final cubes and metadata')
            if split_posang:
                self._logger.debug('> split position angles')
                pas = np.unique(sci_posang)
                for pa in pas:
                    ii = np.where(sci_posang == pa)[0]

                    # save metadata
                    object_files[(object_files['INS4 DROT2 POSANG'] + 90) == pa].to_csv(path.products / 'science_posang={:06.2f}_frames.csv'.format(pa))
                    fits.writeto(path.products / 'science_posang={:06.2f}_posang.fits'.format(pa), sci_posang[ii], overwrite=True)

                    # save final cubes
                    fits.writeto(path.products / 'science_posang={:06.2f}_cube.fits'.format(pa), sci_cube[:, ii], overwrite=True)
            else:
                # save metadata
                object_files.to_csv(path.products / 'science_posang=all_frames.csv')
                fits.writeto(path.products / 'science_posang=all_posang.fits', sci_posang, overwrite=True)

                # save final cubes
                fits.writeto(path.products / 'science_posang=all_cube.fits', sci_cube, overwrite=True)

            # delete big cubes
            del sci_cube

        # update recipe execution
        self._update_recipe_status('sph_ird_combine_data', sphere.SUCCESS)

        # reduction status
        self._status = sphere.COMPLETE


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

        self._logger.info('Clean reduction data')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_clean',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self._path

        # tmp
        if path.tmp.exists():
            self._logger.debug('> remove {}'.format(path.tmp))
            shutil.rmtree(path.tmp, ignore_errors=True)

        # sof
        if path.sof.exists():
            self._logger.debug('> remove {}'.format(path.sof))
            shutil.rmtree(path.sof, ignore_errors=True)

        # calib
        if path.calib.exists():
            self._logger.debug('> remove {}'.format(path.calib))
            shutil.rmtree(path.calib, ignore_errors=True)

        # preproc
        if path.preproc.exists():
            self._logger.debug('> remove {}'.format(path.preproc))
            shutil.rmtree(path.preproc, ignore_errors=True)

        # raw
        if delete_raw:
            if path.raw.exists():
                self._logger.debug('> remove {}'.format(path.raw))
                self._logger.warning('   ==> delete raw files')
                shutil.rmtree(path.raw, ignore_errors=True)

        # products
        if delete_products:
            if path.products.exists():
                self._logger.debug('> remove {}'.format(path.products))
                self._logger.warning('   ==> delete products')
                shutil.rmtree(path.products, ignore_errors=True)

        # update recipe execution
        self._logger.debug('> update recipe execution')
        self._recipes_status['sph_ird_clean'] = True

        # update recipe execution
        self._update_recipe_status('sph_ird_clean', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
