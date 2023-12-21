import pandas as pd
import subprocess
import logging
import numpy as np
import scipy.ndimage as ndimage
import scipy.interpolate as interp
import scipy.optimize as optim
import shutil
import configparser
import collections

from pathlib import Path
from astropy.io import fits
from astropy.modeling import models, fitting

import sphere
import sphere.utils as utils
import sphere.utils.imutils as imutils
import sphere.utils.aperture as aperture
import sphere.utils.toolbox as toolbox
import sphere.utils.transmission as transmission

_log = logging.getLogger(__name__)


class ImagingReduction(object):
    '''
    SPHERE/IRDIS imaging reduction class. It handles both the
    dual-band imaging (DBI) and classical imaging (CI) observing
    modes.
    '''

    ##################################################
    # Class variables
    ##################################################

    # specify for each recipe which other recipes need to have been executed before
    recipe_requirements = collections.OrderedDict([
        ('sort_files', []),
        ('sort_frames', ['sort_files']),
        ('check_files_association', ['sort_files']),
        ('sph_ird_cal_dark', ['sort_files']),
        ('sph_ird_cal_detector_flat', ['sort_files']),
        ('sph_ird_preprocess_science', ['sort_files', 'sort_frames', 'sph_ird_cal_dark',
                                        'sph_ird_cal_detector_flat']),
        ('sph_ird_star_center', ['sort_files', 'sort_frames', 'sph_ird_preprocess_science']),
        ('sph_ird_combine_data', ['sort_files', 'sort_frames', 'sph_ird_preprocess_science']),
        ('sph_ird_clean', [])
    ])

    ##################################################
    # Constructor
    ##################################################

    def __new__(cls, path, clean_start=True, log_level='info', user_config=None, sphere_handler=None):
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

        clean_start : bool
            Remove all results from previous reductions for a clean start.
            Default is True
        
        log_level : {'debug', 'info', 'warning', 'error', 'critical'}
            The log level of the handler

        user_config : str
            Path to a user-provided configuration. Default is None, i.e. the
            reduction will use the package default configuration parameters
        
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
            _log.error(f'No raw/ subdirectory. {path} is not a valid reduction path')
            return None
        else:
            # it's all good: create instance!
            reduction = super(ImagingReduction, cls).__new__(cls)

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

        reduction._logger.info(f'Creating IRDIS imaging reduction at path {path}')

        #
        # v1.4 - True North correction change
        #
        reduction._logger.warning('##################################################################')
        reduction._logger.warning('Since version 1.4 of the pipeline, the default -1.75° true North  ')
        reduction._logger.warning('offset is automatically added to the derotation angles. The offset')
        reduction._logger.warning('value can be modified in the configuration of the reduction:      ')
        reduction._logger.warning('                                                                  ')
        reduction._logger.warning('  >>> reduction.config[\'cal_true_north\'] = xxx                  ')
        reduction._logger.warning('                                                                  ')
        reduction._logger.warning('To avoid any issues, make sure to:                                ')
        reduction._logger.warning('  * either reprocess data previously processed with version <1.4  ')
        reduction._logger.warning('  * or take into account the offset in your astrometric analysis  ')
        reduction._logger.warning('##################################################################')
        
        #
        # clean start
        #
        if clean_start:
            reduction._logger.info('Erase outputs of previous reduction for a clean start')
            reduction._path.remove(delete_raw=False, delete_products=True, logger=reduction._logger)
            config_file = reduction._path.root / 'reduction_config.ini'
            if config_file.exists():
                config_file.unlink()
        
        #
        # configuration
        #
        cfgfile = f'{Path(sphere.__file__).parent}/instruments/{reduction._instrument}.ini'
        cfgparser = configparser.ConfigParser()

        reduction._logger.debug('> read default configuration')
        cfgparser.read(cfgfile)

        # instrument
        reduction._pixel = float(cfgparser.get('instrument', 'pixel'))
        reduction._nwave = 2

        # calibration
        reduction._wave_cal_lasers = np.array(eval(cfgparser.get('calibration', 'wave_cal_lasers')))

        # imaging calibration
        reduction._default_center = np.array(eval(cfgparser.get('calibration-imaging', 'default_center')))
        reduction._orientation_offset = eval(cfgparser.get('calibration-imaging', 'orientation_offset'))

        # reduction parameters
        cfg = {}
        for group in ['reduction', 'reduction-imaging']:
            items = dict(cfgparser.items(group))
            for key, value in items.items():
                try:
                    val = eval(value)
                except NameError:
                    val = value
                cfg[key] = val
        reduction._config = utils.Configuration(reduction._path, reduction._logger, cfg)

        # load user-provided default configuration parameters
        if user_config:
            user_config = Path(user_config).expanduser()

            reduction._config.load_from_file(user_config)
        
        #
        # reduction and recipes status
        #
        reduction._status = sphere.INIT
        reduction._recipes_status = collections.OrderedDict()

        for recipe in reduction.recipe_requirements.keys():
            reduction._update_recipe_status(recipe, sphere.NOTSET)
        
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
        return f'<ImagingReduction, instrument={self._instrument}, mode={self._mode}, path={self._path}, log={self.loglevel}>'

    def __format__(self):
        return self.__repr__()

    ##################################################
    # Properties
    ##################################################

    @property
    def loglevel(self):
        return logging.getLevelName(self._logger.level)

    @loglevel.setter
    def loglevel(self, level):
        self._logger.setLevel(level.upper())
    
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
    def status(self):
        return self._status
    
    @property
    def config(self):
        return self._config

    @property
    def mode(self):
        return self._mode

    ##################################################
    # Generic class methods
    ##################################################

    def init_reduction(self):
        '''
        Sort files and frames, perform sanity check
        '''

        self._logger.info('====> Init <====')

        self.sort_files()
        self.sort_frames()
        self.check_files_association()


    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        self._logger.info('====> Static calibrations <====')

        config = self.config

        self.sph_ird_cal_dark(silent=config['misc_silent_esorex'])
        self.sph_ird_cal_detector_flat(silent=config['misc_silent_esorex'])


    def preprocess_science(self):
        '''
        Clean and collapse images
        '''

        self._logger.info('====> Science pre-processing <====')

        config = self.config

        self.sph_ird_preprocess_science(subtract_background=config['preproc_subtract_background'],
                                        fix_badpix=config['preproc_fix_badpix'],
                                        collapse_science=config['preproc_collapse_science'],
                                        collapse_type=config['preproc_collapse_type'],
                                        coadd_value=config['preproc_coadd_value'],
                                        collapse_psf=config['preproc_collapse_psf'],
                                        collapse_center=config['preproc_collapse_center'])


    def process_science(self):
        '''
        Perform star center, combine cubes into final (x,y,time,lambda)
        cubes, correct anamorphism and scale the images
        '''

        self._logger.info('====> Science processing <====')

        config = self.config

        self.sph_ird_star_center(high_pass_psf=config['center_high_pass_psf'],
                                 high_pass_waffle=config['center_high_pass_waffle'],
                                 offset=config['center_offset'],
                                 box_psf=config['center_box_psf'],
                                 box_waffle=config['center_box_waffle'],
                                 plot=config['misc_plot'])
        self.sph_ird_combine_data(cpix=config['combine_cpix'],
                                  psf_dim=config['combine_psf_dim'],
                                  science_dim=config['combine_science_dim'],
                                  correct_anamorphism=config['combine_correct_anamorphism'],
                                  manual_center=config['combine_manual_center'],
                                  center_selection=config['combine_center_selection'],
                                  coarse_centering=config['combine_coarse_centering'],
                                  shift_method=config['combine_shift_method'],
                                  save_scaled=config['combine_save_scaled'])


    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''

        self._logger.info('====> Clean-up <====')

        config = self.config

        if config['clean']:
            self.sph_ird_clean(delete_raw=config['clean_delete_raw'],
                               delete_products=config['clean_delete_products'],
                               delete_config=config['clean_delete_config'])


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
        path = self.path

        # load existing configuration
        self.config.load()
        
        # files info
        fname = path.preproc / 'files.csv'
        if fname.exists():
            self._logger.debug('> read files.csv')

            files_info = pd.read_csv(fname, index_col=0)

            # convert times
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

            # recipe execution status
            self._update_recipe_status('sort_files', sphere.SUCCESS)
            if np.any(files_info['PRO CATG'] == 'IRD_MASTER_DARK'):
                self._update_recipe_status('sph_ird_cal_dark', sphere.SUCCESS)
            if np.any(files_info['PRO CATG'] == 'IRD_FLAT_FIELD'):
                self._update_recipe_status('sph_ird_cal_detector_flat', sphere.SUCCESS)

            # update instrument mode
            # FutureWarning: Series.__getitem__ treating keys as positions is deprecated.
            # In a future version, integer keys will always be treated as labels
            # (consistent with DataFrame behavior). To access a value by position, use `ser.iloc[pos]`
            # self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]
            self._mode = files_info['INS1 MODE'][files_info['DPR CATG'] == 'SCIENCE'].iloc[0]

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

            # recipe execution status
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

        # additional checks to recipe execution status
        if frames_info_preproc is not None:
            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = f'{file}_DIT{idx:03d}_preproc'
                file = list(path.preproc.glob(f'{fname}.fits'))
                done = done and (len(file) == 1)
            if done:
                self._update_recipe_status('sph_ird_preprocess_science', sphere.SUCCESS)
            self._logger.debug(f'> sph_ird_preprocess_science status = {done}')

            done = True
            files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
                                        (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
            for file, idx in files:
                fname = f'{file}_DIT{idx:03d}_preproc_centers'
                file = list(path.preproc.glob(f'{fname}.fits'))
                done = done and (len(file) == 1)
            if done:
                self._update_recipe_status('sph_ird_star_center', sphere.SUCCESS)
            self._logger.debug(f'> sph_ird_star_center status = {done}')
            
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

        # recipe execution status
        self._update_recipe_status('sort_files', sphere.NOTSET)

        # parameters
        path = self.path

        # list files
        files = path.raw.glob('*.fits')
        files = [f.stem for f in files]

        if len(files) == 0:
            self._logger.critical('No raw FITS files in reduction path')
            self._update_recipe_status('sort_files', sphere.ERROR)
            self._status = sphere.FATAL
            return

        self._logger.info(f' * found {len(files)} raw FITS files')

        # read list of keywords
        self._logger.debug('> read keyword list')
        keywords = []
        file = open(Path(sphere.__file__).parent / 'instruments' / 'keywords_irdifs.dat', 'r')
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
        # set dtype to object because there is a mix of float and str columns
        # any float columns will be inferred from the data (see infer_objects)
        files_info = pd.DataFrame(index=pd.Index(files, name='FILE'), columns=keywords_short, dtype=object)

        self._logger.debug('> read FITS keywords')
        for f in files:
            hdu = fits.open(path.raw / f'{f}.fits')
            hdr = hdu[0].header

            for k, sk in zip(keywords, keywords_short):
                if k == 'HIERARCH ESO INS4 DROT2 BEGIN':
                    # in June 2021 ESO changed INS4 DROT2 BEGIN to INS4 DROT2 START
                    v_begin = hdr.get('HIERARCH ESO INS4 DROT2 BEGIN')
                    v_start = hdr.get('HIERARCH ESO INS4 DROT2 START')
                    files_info.loc[f, sk] = v_begin if v_begin else v_start
                elif k == 'HIERARCH ESO INS4 DROT3 BEGIN':
                    # in June 2021 ESO changed INS4 DROT3 BEGIN to INS4 DROT3 START
                    v_begin = hdr.get('HIERARCH ESO INS4 DROT3 BEGIN')
                    v_start = hdr.get('HIERARCH ESO INS4 DROT3 START')
                    files_info.loc[f, sk] = v_begin if v_begin else v_start
                else:
                    files_info.loc[f, sk] = hdr.get(k)

            hdu.close()

        # drop files that are not handled, based on DPR keywords
        self._logger.debug('> drop unsupported file types')
        files_info.dropna(subset=['DPR TYPE'], inplace=True)
        files_info = files_info[(files_info['DPR CATG'] != 'ACQUISITION') & (files_info['DPR TYPE'] != 'OBJECT,AO')]

        # check instruments
        instru = files_info['SEQ ARM'].unique()
        if len(instru) != 1:
            self._logger.critical(f'Sequence is mixing different instruments: {instru}')
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
        # FutureWarning: Series.__getitem__ treating keys as positions is deprecated.
        # In a future version, integer keys will always be treated as labels
        # (consistent with DataFrame behavior). To access a value by position, use `ser.iloc[pos]`
        # self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]
        self._mode = files_info['INS1 MODE'][files_info['DPR CATG'] == 'SCIENCE'].iloc[0]

        # determine data types for columns with float instead of object
        files_info = files_info.infer_objects()

        # sort by acquisition time
        files_info.sort_values(by='DATE-OBS', inplace=True)

        # save files_info
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')
        self._files_info = files_info

        # recipe execution status
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
        path = self.path
        files_info = self.files_info

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
        true_north = self.config['cal_true_north']        
        ret = toolbox.compute_angles(frames_info, true_north, logger=self._logger)
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

        ra_drot   = cinfo['INS4 DROT2 RA'].iloc[0]
        ra_drot_h = np.floor(ra_drot/1e4)
        ra_drot_m = np.floor((ra_drot - ra_drot_h*1e4)/1e2)
        ra_drot_s = ra_drot - ra_drot_h*1e4 - ra_drot_m*1e2
        RA = f'{ra_drot_h:02.0f}:{ra_drot_m:02.0f}:{ra_drot_s:02.3f}'

        dec_drot  = cinfo['INS4 DROT2 DEC'].iloc[0]
        sign = np.sign(dec_drot)
        udec_drot  = np.abs(dec_drot)
        dec_drot_d = np.floor(udec_drot/1e4)
        dec_drot_m = np.floor((udec_drot - dec_drot_d*1e4)/1e2)
        dec_drot_s = udec_drot - dec_drot_d*1e4 - dec_drot_m*1e2
        dec_drot_d *= sign
        DEC = f'{dec_drot_d:02.0f}:{dec_drot_m:02.0f}:{dec_drot_s:02.2f}'

        pa_start = cinfo['PARANG'].iloc[0]
        pa_end   = cinfo['PARANG'].iloc[-1]

        posang  = cinfo['INS4 DROT2 POSANG'].unique()
        posangs = [f'{p:.2f}°' for p in posang]
        
        date = str(cinfo['DATE'].iloc[0])[0:10]
        
        self._logger.info(f" * Programme ID: {cinfo['OBS PROG ID'].iloc[0]}")
        self._logger.info(f" * OB name:      {cinfo['OBS NAME'].iloc[0]}")
        self._logger.info(f" * OB ID:        {cinfo['OBS ID'].iloc[0]}")
        self._logger.info(f" * Object:       {cinfo['OBJECT'].iloc[0]}")
        self._logger.info(f' * RA / DEC:     {RA} / {DEC}')
        self._logger.info(f' * Date:         {date}')
        self._logger.info(f" * Instrument:   {cinfo['SEQ ARM'].iloc[0]}")
        self._logger.info(f" * Derotator:    {cinfo['INS4 DROT2 MODE'].iloc[0]}")
        self._logger.info(f" * VIS WFS mode: {cinfo['AOS VISWFS MODE'].iloc[0]}")
        self._logger.info(f" * IR WFS mode:  {cinfo['AOS IRWFS MODE'].iloc[0]}")
        self._logger.info(f" * Coronagraph:  {cinfo['INS COMB ICOR'].iloc[0]}")
        self._logger.info(f" * Mode:         {cinfo['INS1 MODE'].iloc[0]}")
        self._logger.info(f" * Filter:       {cinfo['INS COMB IFLT'].iloc[0]}")
        self._logger.info(f" * DIT:          {cinfo['DET SEQ1 DIT'].iloc[0]:.2f} sec")
        self._logger.info(f" * NDIT:         {cinfo['DET NDIT'].iloc[0]:.0f}")
        self._logger.info(f" * Texp:         {cinfo['DET SEQ1 DIT'].sum() / 60:.2f} min")
        self._logger.info(f' * PA:           {pa_start:.2f}° ==> {pa_end:.2f}° = {np.abs(pa_end - pa_start):.2f}°')
        self._logger.info(f" * POSANG:       {', '.join(posangs)}")

        # recipe execution status
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
        files_info = self.files_info

        # instrument arm
        arm = files_info['SEQ ARM'].unique()
        if len(arm) != 1:
            self._logger.error(f'Sequence is mixing different instruments: {arm}')
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        # IRDIS obs mode and filter combination
        modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'].unique()
        if len(modes) != 1:
            self._logger.error(f'Sequence is mixing different types of observations: {modes}')
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return

        filter_combs = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()
        if len(filter_combs) != 1:
            self._logger.error(f'Sequence is mixing different types of filters combinations: {filter_combs}')
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return
        filter_comb = filter_combs[0]

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
            self._logger.error(f' * there should be more than 1 flat in filter combination {filter_comb}')

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
                self._logger.warning(f' * there is no dark/background for science files with DIT={DIT} sec. It is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive')

            # sky backgrounds
            cfiles = files_info[(files_info['DPR TYPE'] == 'SKY') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
            if len(cfiles) == 0:
                warning_flag += 1
                self._logger.warning(f' * there is no sky background for science files with DIT={DIT} sec. Using a sky background instead of an internal instrumental background can usually provide a cleaner data reduction, especially in K-band')

        # error reporting
        self._logger.debug('> report status')
        if error_flag:
            self._logger.error(f'There are {warning_flag} warning(s) and {error_flag} error(s) in the classification of files')
            self._update_recipe_status('check_files_association', sphere.ERROR)
            return
        else:
            self._logger.warning(f'There are {warning_flag} warning(s) and {error_flag} error(s) in the classification of files')

        # recipe execution status
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
        path = self.path
        files_info = self.files_info

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

                    self._logger.info(f' * {ctype} in filter {cfilt} with DIT={DIT:.2f} sec ({len(cfiles)} files)')

                    # create sof
                    self._logger.debug('> create sof file')
                    sof = path.sof / f'dark_filt={cfilt}_DIT={DIT:.2f}.sof'
                    file = open(sof, 'w')
                    for f in files:
                        file.write(f"{path.raw}/{f}.fits     IRD_DARK_RAW\n")
                    file.close()

                    # products
                    if ctype == 'SKY':
                        loc = 'sky'
                    else:
                        loc = 'internal'
                    dark_file = f'dark_{loc}_filt={cfilt}_DIT={DIT:.2f}'
                    bpm_file  = f'dark_{loc}_bpm_filt={cfilt}_DIT={DIT:.2f}'

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
                            f'--ird.master_dark.max_acceptable={max_level}',
                            f'--ird.master_dark.outfilename={path.calib}/{dark_file}.fits',
                            f'--ird.master_dark.badpixfilename={path.calib}/{bpm_file}.fits',
                            str(sof)]

                    # check esorex
                    if shutil.which('esorex') is None:
                        self._logger.error('esorex does not appear to be in your PATH. Please make sure that the ESO pipeline is properly installed before running vlt-sphere.')
                        self._update_recipe_status('sph_ird_cal_dark', sphere.ERROR)
                        return

                    # execute esorex
                    self._logger.debug(f"> execute {' '.join(args)}")
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
                    files_info.loc[dark_file, 'DPR CATG'] = cfiles['DPR CATG'].iloc[0]
                    files_info.loc[dark_file, 'DPR TYPE'] = cfiles['DPR TYPE'].iloc[0]
                    files_info.loc[dark_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'].iloc[0]
                    files_info.loc[dark_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'].iloc[0]
                    files_info.loc[dark_file, 'INS1 MODE'] = cfiles['INS1 MODE'].iloc[0]
                    files_info.loc[dark_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'].iloc[0]
                    files_info.loc[dark_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'].iloc[0]
                    files_info.loc[dark_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'].iloc[0]
                    files_info.loc[dark_file, 'PROCESSED'] = True
                    files_info.loc[dark_file, 'PRO CATG'] = 'IRD_MASTER_DARK'

                    files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'].iloc[0]
                    files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'].iloc[0]
                    files_info.loc[bpm_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'].iloc[0]
                    files_info.loc[bpm_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'].iloc[0]
                    files_info.loc[bpm_file, 'INS1 MODE'] = cfiles['INS1 MODE'].iloc[0]
                    files_info.loc[bpm_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'].iloc[0]
                    files_info.loc[bpm_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'].iloc[0]
                    files_info.loc[bpm_file, 'PROCESSED'] = True
                    files_info.loc[bpm_file, 'PRO CATG']  = 'IRD_STATIC_BADPIXELMAP'

        # save
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')

        # recipe execution status
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
        path = self.path
        files_info = self.files_info

        # get list of files
        calibs = files_info[np.logical_not(files_info['PROCESSED']) &
                            ((files_info['DPR TYPE'] == 'FLAT,LAMP') |
                             (files_info['DPR TECH'] == 'IMAGE'))]
        filter_combs = calibs['INS COMB IFLT'].unique()

        for cfilt in filter_combs:
            cfiles = calibs[calibs['INS COMB IFLT'] == cfilt]
            files = cfiles.index

            self._logger.info(f' * filter {cfilt} ({len(cfiles)} files)')

            # create sof
            self._logger.debug('> create sof file')
            sof = path.sof / f'flat_filt={cfilt}.sof'
            file = open(sof, 'w')
            for f in files:
                file.write(f"{path.raw}/{f}.fits     IRD_FLAT_FIELD_RAW\n")
            file.close()

            # products
            flat_file = f'flat_filt={cfilt}'
            bpm_file  = f'flat_bpm_filt={cfilt}'

            # esorex parameters
            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_instrument_flat',
                    '--ird.instrument_flat.save_addprod=TRUE',
                    f'--ird.instrument_flat.outfilename={path.calib}/{flat_file}.fits',
                    f'--ird.instrument_flat.badpixfilename={path.calib}/{bpm_file}.fits',
                    str(sof)]

            # check esorex
            if shutil.which('esorex') is None:
                self._logger.error('esorex does not appear to be in your PATH. Please make sure that the ESO pipeline is properly installed before running vlt-sphere.')
                self._update_recipe_status('sph_ird_cal_detector_flat', sphere.ERROR)
                return

            # execute esorex
            self._logger.debug(f"> execute {' '.join(args)}")
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
            files_info.loc[flat_file, 'DPR CATG'] = cfiles['DPR CATG'].iloc[0]
            files_info.loc[flat_file, 'DPR TYPE'] = cfiles['DPR TYPE'].iloc[0]
            files_info.loc[flat_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'].iloc[0]
            files_info.loc[flat_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'].iloc[0]
            files_info.loc[flat_file, 'INS1 MODE'] = cfiles['INS1 MODE'].iloc[0]
            files_info.loc[flat_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'].iloc[0]
            files_info.loc[flat_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'].iloc[0]
            files_info.loc[flat_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'].iloc[0]
            files_info.loc[flat_file, 'PROCESSED'] = True
            files_info.loc[flat_file, 'PRO CATG'] = 'IRD_FLAT_FIELD'

            files_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'].iloc[0]
            files_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'].iloc[0]
            files_info.loc[bpm_file, 'INS COMB IFLT'] = cfiles['INS COMB IFLT'].iloc[0]
            files_info.loc[bpm_file, 'INS4 FILT2 NAME'] = cfiles['INS4 FILT2 NAME'].iloc[0]
            files_info.loc[bpm_file, 'INS1 MODE'] = cfiles['INS1 MODE'].iloc[0]
            files_info.loc[bpm_file, 'INS1 FILT NAME'] = cfiles['INS1 FILT NAME'].iloc[0]
            files_info.loc[bpm_file, 'INS1 OPTI2 NAME'] = cfiles['INS1 OPTI2 NAME'].iloc[0]
            files_info.loc[bpm_file, 'PROCESSED'] = True
            files_info.loc[bpm_file, 'PRO CATG']  = 'IRD_NON_LINEAR_BADPIXELMAP'

        # save
        self._logger.debug('> save files.csv')
        files_info.to_csv(path.preproc / 'files.csv')

        # recipe execution status
        self._update_recipe_status('sph_ird_cal_detector_flat', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_preprocess_science(self,
                                   subtract_background=True, fix_badpix=True,
                                   collapse_science=False, collapse_type='mean', coadd_value=2,
                                   collapse_psf=True, collapse_center=True):
        '''Pre-processes the science frames.

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

        self._logger.info('Pre-process science files')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_preprocess_science',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info = self.files_info
        frames_info = self.frames_info

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
            bpm_files = [path.calib / f'{f}.fits' for f in bpm_files]
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
            self._logger.error(f'There should be exactly 1 flat file. Found {len(flat_file)}.')
            self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
            return
        flat = fits.getdata(path.calib / f'{flat_file.index[0]}.fits')

        # final dataframe
        self._logger.debug('> create frames_info_preproc data frame')
        index = pd.MultiIndex(names=['FILE', 'IMG'], levels=[[], []], codes=[[], []])
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

                self._logger.info(f'{len(sfiles)} files of type {typ} with DIT={DIT} sec')

                if subtract_background:
                    # look for sky, then background, then darks
                    # normally there should be only one with the proper DIT
                    dfiles = []
                    for d in dark_types:
                        dfiles = files_info[(files_info['PRO CATG'] == 'IRD_MASTER_DARK') &
                                            (files_info['DPR TYPE'] == d) & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
                        if len(dfiles) != 0:
                            break
                    self._logger.info(f'   ==> found {len(dfiles)} corresponding {d} file')

                    if len(dfiles) == 0:
                        # issue a warning if absolutely no background is found
                        self._logger.warning('No background has been found. Pre-processing will continue but data quality will likely be affected')
                        bkg = np.zeros((1024, 2048))
                    elif len(dfiles) == 1:
                        bkg = fits.getdata(path.calib / f'{dfiles.index[0]}.fits')
                    elif len(dfiles) > 1:
                        # FIXME: handle cases when multiple backgrounds are found?
                        self._logger.error(f'Unexpected number of background files ({len(dfiles)})')
                        self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                        return

                # process files
                for idx, (fname, finfo) in enumerate(sfiles.iterrows()):
                    # frames_info extract
                    finfo = frames_info.loc[(fname, slice(None)), :]

                    self._logger.info(f' * file {idx + 1}/{len(sfiles)}: {fname}, NDIT={len(finfo)}')

                    # read data
                    self._logger.info('   ==> read data')
                    img, hdr = fits.getdata(path.raw / f'{fname}.fits', header=True)

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
                    true_north = self.config['cal_true_north']
                    if (typ == 'OBJECT,CENTER'):
                        if collapse_center:
                            self._logger.info('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'mean', logger=self._logger)
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'none', logger=self._logger)
                    elif (typ == 'OBJECT,FLUX'):
                        if collapse_psf:
                            self._logger.info('   ==> collapse: mean')
                            img = np.mean(img, axis=0, keepdims=True)
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'mean', logger=self._logger)
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'none', logger=self._logger)
                    elif (typ == 'OBJECT'):
                        if collapse_science:
                            if collapse_type == 'mean':
                                self._logger.info(f'   ==> collapse: mean ({len(img)} -> 1 frame, 0 dropped)')
                                img = np.mean(img, axis=0, keepdims=True)

                                frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'mean', logger=self._logger)
                            elif collapse_type == 'coadd':
                                if (not isinstance(coadd_value, int)) or (coadd_value <= 1):
                                    self._logger.error('coadd_value must be an integer >1')
                                    self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                                    return

                                coadd_value = int(coadd_value)
                                NDIT = len(img)
                                NDIT_new = NDIT // coadd_value
                                dropped = NDIT % coadd_value

                                if coadd_value > NDIT:
                                    self._logger.error(f'coadd_value ({coadd_value}) must be < NDIT ({NDIT})')
                                    self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                                    return

                                self._logger.info(f'   ==> collapse: coadd by {coadd_value} ({NDIT} -> {NDIT_new} frames, {dropped} dropped)')

                                # coadd frames
                                nimg = np.empty((NDIT_new, 1024, 2048), dtype=img.dtype)
                                for f in range(NDIT_new):
                                    nimg[f] = np.mean(img[f*coadd_value:(f+1)*coadd_value], axis=0)
                                img = nimg

                                frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'coadd', coadd_value=coadd_value, logger=self._logger)
                            else:
                                self._logger.error(f'Unknown collapse type {collapse_type}')
                                self._update_recipe_status('sph_ird_preprocess_science', sphere.ERROR)
                                return
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, true_north, 'none', logger=self._logger)

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
                        fits.writeto(path.preproc / f'{fname}_DIT{f:03d}_preproc.fits', frame, hdr,
                                     overwrite=True, output_verify='silentfix')

        # sort and save final dataframe
        self._logger.debug('> save frames_info_preproc.csv')
        frames_info_preproc.sort_values(by='TIME', inplace=True)
        frames_info_preproc.to_csv(path.preproc / 'frames_preproc.csv')

        self._frames_info_preproc = frames_info_preproc

        # recipe execution status
        self._update_recipe_status('sph_ird_preprocess_science', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_star_center(self, high_pass_psf=False, high_pass_waffle=False, offset=(0, 0), box_psf=60, box_waffle=16, plot=True):
        '''Determines the star center for all frames where a center can be
        determined (OBJECT,CENTER and OBJECT,FLUX)

        Parameters
        ----------
        high_pass_psf : bool
            Apply high-pass filter to the PSF image before searching for the center.
            Default is False

        high_pass_waffle : bool
            Apply high-pass filter to the center image before searching for the waffle spots.
            Default is False

        offset : tuple
            Apply an (x,y) offset to the default center position, for the waffle centering.
            The offset will move the search box of the waffle spots by the amount of
            specified pixels in each direction. Default is no offset

        box_psf : int
            Size of the box in which the PSF fit is performed. Default is 60 pixels

        box_waffle : int
            Size of the box in which the waffle fit is performed. Default is 16 pixels

        plot : bool
            Display and save diagnostic plot for quality check. Default is True

        '''

        self._logger.info('Star centers determination')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_star_center',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        pixel = self.pixel
        orientation_offset = self._orientation_offset
        center_guess = self._default_center
        frames_info = self.frames_info_preproc

        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)
        wave = np.array(wave)

        # start with OBJECT,FLUX
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        if len(flux_files) != 0:
            for file, idx in flux_files.index:
                self._logger.info(f' * OBJECT,FLUX: {file}')

                # read data
                self._logger.debug('> read data')
                fname = f'{file}_DIT{idx:03d}_preproc'
                cube, hdr = fits.getdata(path.preproc / f'{fname}.fits', header=True)

                # centers
                if plot:
                    save_path = path.products / f'{fname}_psf_fitting.pdf'
                else:
                    save_path = None
                img_center = toolbox.star_centers_from_PSF_img_cube(cube, wave, pixel, exclude_fraction=0.3,
                                                                    high_pass=high_pass_psf, box_size=box_psf,
                                                                    save_path=save_path, logger=self._logger)

                # save
                self._logger.debug('> save centers')
                fits.writeto(path.preproc / f'{fname}_centers.fits', img_center, overwrite=True)

        # then OBJECT,CENTER
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) != 0:
            for file, idx in starcen_files.index:
                self._logger.info(f' * OBJECT,CENTER: {file}')

                # read data
                self._logger.debug('> read data')
                fname = f'{file}_DIT{idx:03d}_preproc'
                cube, hdr = fits.getdata(path.preproc / f'{fname}.fits', header=True)

                # coronagraph
                coro_name = starcen_files.loc[(file, idx), 'INS COMB ICOR']
                if coro_name == 'N_NS_CLEAR':
                    coro = False
                else:
                    coro = True

                # centers
                waffle_orientation = hdr['HIERARCH ESO OCS WAFFLE ORIENT']
                self._logger.debug(f'> waffle orientation: {waffle_orientation}')
                if plot:
                    save_path = path.products / f'{fname}_waffle_fitting.pdf'
                else:
                    save_path = None
                spot_center, spot_dist, img_center \
                    = toolbox.star_centers_from_waffle_img_cube(cube, wave, waffle_orientation, center_guess,
                                                                pixel, orientation_offset, high_pass=high_pass_waffle,
                                                                center_offset=offset, box_size=box_waffle,
                                                                coro=coro, save_path=save_path, logger=self._logger)

                # save
                self._logger.debug('> save centers')
                fits.writeto(path.preproc / f'{fname}_centers.fits', img_center, overwrite=True)

        # recipe execution status
        self._update_recipe_status('sph_ird_star_center', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE


    def sph_ird_combine_data(self, cpix=True, psf_dim=80, science_dim=290, correct_anamorphism=True,
                             shift_method='fft', manual_center=None, center_selection='first',
                             coarse_centering=False, save_scaled=False):
        '''Combine and save the science data into final cubes

        All types of data are combined independently: PSFs
        (OBJECT,FLUX), star centers (OBJECT,CENTER) and standard
        coronagraphic images (OBJECT). For each type of data, the
        method saves 4 or 5 different files:

          - *_cube: the (x,y,time,lambda) cube

          - *_derot: the derotation angles vector. This vector takes
                     into account the parallactic angle, the default 
                     -1.75° true North offset, and any instrumental 
                     pupil offset. This is the values that need to be 
                     used for aligning the images with North up and 
                     East left.

          - *_frames: a csv file with all the information for every
                      frames. There is one line by time step in the
                      data cube.

          - *_cube_scaled: the (x,y,time,lambda) cube with images
                           rescaled spectraly. This is useful if you
                           plan to perform spectral differential
                           imaging in your analysis.

        Centering
        ---------

        By default, a fine (sub-pixel) centering is performed if the
        an OBJECT,CENTER frame was acquired in the sequence or if
        there is a valid user-provided center. However, if the
        coarse_centering keyword is set to True, only a "coarse
        centering" is performed, which requires no interpolation:

          - only integer shifts (shift_method='roll')
          - centering on an integer pixel (cpix=True)
          - no correction of the anamorphism (correct_anamorphism=False)
          - no saving of the rescaled frames (save_scaled=False)

        This option is useful if the user wants to perform a
        posteriori centering of the frames, e.g. to fully preserve
        photometry.

        If there was no OBJECT,CENTER acquired in the sequence, then
        the centering will be performed with respect to a default,
        pre-defined center that is representative of the typical
        center of the coronagraph.

        Parameters
        ----------
        cpix : bool
            If True the images are centered on the pixel at coordinate
            (dim//2,dim//2). If False the images are centered between
            4 pixels, at coordinates ((dim-1)/2,(dim-1)/2). The value
            of cpix is automatically set to True when coarse_centering
            is set to True. Default is True.

        psf_dim : even int
            Size of the PSF images. Default is 80x80 pixels

        science_dim : even int
            Size of the science images (star centers and standard
            coronagraphic images). Default is 290, 290 pixels

        correct_anamorphism : bool
            Correct the optical anamorphism of the instrument (see
            user manual for details). The value of correct_anamorphism
            is automatically set to True when coarse_centering is set
            to True. Default is True.

        manual_center : array
            User provided centers for the OBJECT,CENTER and OBJECT
            frames. This should be an array of either 2 or nwave*2
            values. Default is None

        center_selection : str        
            Specify which star center to use when multiple are
            available. Possible values are first, last, and time. The
            time option indicates to use the star center file that is
            closest in time with respect to each science file. Default
            is first
        
        coarse_centering : bool
            Control if images are finely centered or not before being
            combined. However the images are still roughly centered by
            shifting them by an integer number of pixel to bring the
            center of the data close to the center of the images. This
            option is useful if fine centering must be done
            afterwards. Default is False.

        shift_method : str
            Method to scaling and shifting the images: fft or interp.
            Default is fft

        save_scaled : bool
            Also save the wavelength-rescaled cubes. Makes the process
            much longer. The value of save_scaled is automatically set
            to False when coarse_centering is set to True. The default
            is False
        '''

        self._logger.info('Combine science data')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_combine_data',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        nwave = self.nwave
        frames_info = self.frames_info_preproc

        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)
        wave = np.array(wave)

        self._logger.debug('> save final wavelength')
        fits.writeto(path.products / 'wavelength.fits', wave, overwrite=True)

        # max images size
        if psf_dim > 1024:
            self._logger.warning('psf_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            psf_dim = 1024

        if science_dim > 1024:
            self._logger.warning('science_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            science_dim = 1024

        # centering configuration
        if coarse_centering:
            self._logger.warning('Images will be coarsely centered without any interpolation. Automatic settings for coarse centering: shift_method=\'roll\', cpix=True, correct_anamorphism=False, save_scaled=False')
            shift_method = 'roll'
            cpix = True
            correct_anamorphism = False
            save_scaled = False

        if manual_center is not None:
            manual_center = np.array(manual_center)

            if (manual_center.shape != (2,)) and (manual_center.shape != (nwave, 2)):
                self._logger.error('manual_center does not have the right number of dimensions.')
                self._update_recipe_status('sph_ird_combine_data', sphere.ERROR)
                return

            if manual_center.shape == (2,):
                manual_center = np.full((nwave, 2), manual_center, dtype=np.float64)

            self._logger.warning('Images will be centered using the user-provided center ({},{})'.format(*manual_center[0]))

        #
        # OBJECT,FLUX
        #
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        nfiles = len(flux_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT,FLUX data')

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
                self._logger.info(f'   ==> file {file_idx + 1}/{len(flux_files)}: {file}, DIT #{idx}')

                # read data
                self._logger.debug('> read data')
                fname = f'{file}_DIT{idx:03d}_preproc'
                cube = fits.getdata(path.preproc / f'{fname}.fits')

                self._logger.debug('> read centers')
                cfile = path.preproc / f'{fname}_centers.fits'
                if cfile.exists():
                    centers = fits.getdata(cfile)
                else:
                    self._logger.warning('sph_ird_star_center() has not been executed. Images will be centered using default center ({},{})'.format(*self._default_center))
                    centers = self._default_center

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)

                # neutral density
                self._logger.debug('> read neutral density information')
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                psf_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                psf_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    self._logger.debug(f'> wave {wave_idx}')
                    cx, cy = centers[wave_idx, :]

                    self._logger.debug('> shift and normalize')
                    img  = img.astype(np.float64)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    psf_cube[wave_idx, file_idx] = nimg[:psf_dim, :psf_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        self._logger.debug('> correct anamorphism')
                        nimg = psf_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        psf_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        self._logger.debug('> spatial scaling')
                        nimg = psf_cube[wave_idx, file_idx]
                        psf_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            self._logger.debug('> save final cubes and metadata')
            flux_files.to_csv(path.products / 'psf_frames.csv')
            fits.writeto(path.products / 'psf_cube.fits', psf_cube, overwrite=True)
            fits.writeto(path.products / 'psf_derot.fits', psf_derot, overwrite=True)
            if save_scaled:
                self._logger.debug('> save scaled cubes')
                fits.writeto(path.products / 'psf_cube_scaled.fits', psf_cube_scaled, overwrite=True)

            # delete big cubes
            self._logger.debug('> free memory')
            del psf_cube
            if save_scaled:
                del psf_cube_scaled

        #
        # OBJECT,CENTER
        #
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        nfiles = len(starcen_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT,CENTER data')

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
                self._logger.info(f'   ==> file {file_idx + 1}/{len(starcen_files)}: {file}, DIT #{idx}')

                # read data
                self._logger.debug('> read data')
                fname = f'{file}_DIT{idx:03d}_preproc'
                cube = fits.getdata(path.preproc / f'{fname}.fits')

                # use manual center if explicitely requested
                self._logger.debug('> read centers')
                if manual_center is not None:
                    centers = manual_center
                else:
                    # otherwise read center data
                    centers = fits.getdata(path.preproc / f'{fname}_centers.fits')

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)

                # neutral density
                self._logger.debug('> read neutral density information')
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                cen_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                cen_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    self._logger.debug(f'> wave {wave_idx}')
                    cx, cy = centers[wave_idx, :]

                    self._logger.debug('> shift and normalize')
                    img  = img.astype(np.float64)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    cen_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        self._logger.debug('> correct anamorphism')
                        nimg = cen_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        cen_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        self._logger.debug('> spatial scaling')
                        nimg = cen_cube[wave_idx, file_idx]
                        cen_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            self._logger.debug('> save final cubes and metadata')
            starcen_files.to_csv(path.products / 'starcenter_frames.csv')
            fits.writeto(path.products / 'starcenter_cube.fits', cen_cube, overwrite=True)
            fits.writeto(path.products / 'starcenter_derot.fits', cen_derot, overwrite=True)
            if save_scaled:
                self._logger.debug('> save scaled cubes')
                fits.writeto(path.products / 'starcenter_cube_scaled.fits', cen_cube_scaled, overwrite=True)

            # delete big cubes
            self._logger.debug('> free memory')
            del cen_cube
            if save_scaled:
                del cen_cube_scaled

        #
        # OBJECT
        #
        object_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT']
        nfiles = len(object_files)
        if nfiles != 0:
            self._logger.info(' * OBJECT data')

            # null value for Dithering Motion Stage by default
            dms_dx_ref = 0
            dms_dy_ref = 0

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
                self._logger.info(f'   ==> file {file_idx + 1}/{len(object_files)}: {file}, DIT #{idx}')

                # use manual center if explicitely requested
                self._logger.debug('> read centers')
                if manual_center is not None:
                    centers = manual_center
                else:
                    # otherwise, look whether we have an OBJECT,CENTER frame and select the one requested by user
                    starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
                    if len(starcen_files) == 0:
                        self._logger.warning('No OBJECT,CENTER file in the dataset. Images will be centered using default center ({},{})'.format(*self._default_center))
                        centers = self._default_center
                    else:
                        # selection of the proper OBJECT,CENTER
                        center_selection = center_selection.lower()
                        if center_selection == 'first':
                            center_index = 0
                        elif center_selection == 'last':
                            center_index = len(starcen_files.index.values)-1
                        elif center_selection == 'time':
                            time_cen = starcen_files['DATE-OBS']
                            time_sci = frames_info.loc[(file, idx), 'DATE-OBS']
                            center_index = np.abs(time_sci - time_cen).argmin()
                        else:
                            self._logger.error(f'Unknown OBJECT,CENTER selection {center_selection}. Possible values are first, last, and time.')
                            self._update_recipe_status('sph_ird_combine_data', sphere.ERROR)
                            return

                        fname = f'{starcen_files.index.values[center_index][0]}_DIT{starcen_files.index.values[center_index][1]:03d}_preproc_centers.fits'
                        fpath = path.preproc / fname
                        if fpath.exists():
                            centers = fits.getdata(fpath)

                            # Dithering Motion Stage for star center: value is in micron,
                            # and the pixel size is 18 micron
                            dms_dx_ref = starcen_files['INS1 PAC X'].iloc[0] / 18
                            dms_dy_ref = starcen_files['INS1 PAC Y'].iloc[0] / 18
                        else:
                            self._logger.warning('sph_ird_star_center() has not been executed. Images will be centered using default center ({},{})'.format(*self._default_center))
                            centers = self._default_center

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)
                    dms_dx_ref = np.int(dms_dx_ref)
                    dms_dy_ref = np.int(dms_dy_ref)

                # read data
                self._logger.debug('> read data')
                fname = f'{file}_DIT{idx:03d}_preproc'
                files = list(path.preproc.glob(f'{fname}*.fits'))
                cube = fits.getdata(files[0])

                # neutral density
                self._logger.debug('> read neutral density information')
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                self._logger.debug('> read angles')
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                sci_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                sci_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # Dithering Motion Stage for star center: value is in micron,
                # and the pixel size is 18 micron
                self._logger.debug('> read DMS position')
                dms_dx = frames_info.loc[(file, idx), 'INS1 PAC X'] / 18
                dms_dy = frames_info.loc[(file, idx), 'INS1 PAC Y'] / 18

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    dms_dx = np.int(dms_dx)
                    dms_dy = np.int(dms_dy)

                # center frames
                for wave_idx, img in enumerate(cube):
                    self._logger.debug(f'> wave {wave_idx}')
                    cx, cy = centers[wave_idx, :]

                    # DMS contribution
                    cx = cx + dms_dx_ref + dms_dx
                    cy = cy + dms_dy_ref + dms_dy

                    self._logger.debug('> shift and normalize')
                    img  = img.astype(np.float64)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    sci_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        self._logger.debug('> correct anamorphism')
                        nimg = sci_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        sci_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        self._logger.debug('> spatial scaling')
                        nimg = sci_cube[wave_idx, file_idx]
                        sci_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            self._logger.debug('> save final cubes and metadata')
            object_files.to_csv(path.products / 'science_frames.csv')
            fits.writeto(path.products / 'science_cube.fits', sci_cube, overwrite=True)
            fits.writeto(path.products / 'science_derot.fits', sci_derot, overwrite=True)
            if save_scaled:
                self._logger.debug('> save scaled cubes')
                fits.writeto(path.products / 'science_cube_scaled.fits', sci_cube_scaled, overwrite=True)

            # delete big cubes
            self._logger.debug('> free memory')
            del sci_cube
            if save_scaled:
                del sci_cube_scaled

        # recipe execution status
        self._update_recipe_status('sph_ird_combine_data', sphere.SUCCESS)

        # reduction status
        self._status = sphere.COMPLETE


    def sph_ird_clean(self, delete_raw=False, delete_products=False, delete_config=False):
        '''
        Clean everything except for raw data and science products (by default)

        Parameters
        ----------
        delete_raw : bool
            Delete raw data. Default is False

        delete_products : bool
            Delete science products. Default is False

        delete_config : bool
            Delete configuration file. Default is False
        '''

        self._logger.info('Clean reduction data')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_ird_clean',
                                         self.recipe_requirements, logger=self._logger):
            return

        # remove sub-directories
        self.path.remove(delete_raw=delete_raw, delete_products=delete_products, logger=self._logger)

        # remove config
        if delete_config:
            self.config._file.unlink()
        
        # recipe execution status
        self._update_recipe_status('sph_ird_clean', sphere.SUCCESS)

        # reduction status
        self._status = sphere.COMPLETE
        
