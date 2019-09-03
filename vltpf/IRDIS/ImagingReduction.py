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

from pathlib import Path
from astropy.io import fits
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages

import vltpf
import vltpf.utils as utils
import vltpf.utils.imutils as imutils
import vltpf.utils.aperture as aperture
import vltpf.transmission as transmission
import vltpf.toolbox as toolbox


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
    recipe_requirements = {
        'sort_files': [],
        'sort_frames': ['sort_files'],
        'check_files_association': ['sort_files'],
        'sph_ird_cal_dark': ['sort_files'],
        'sph_ird_cal_detector_flat': ['sort_files'],
        'sph_ird_preprocess_science': ['sort_files', 'sort_frames', 'sph_ird_cal_dark', 
                                       'sph_ird_cal_detector_flat'],
        'sph_ird_star_center': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science'],
        'sph_ird_combine_data': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science']
    }

    ##################################################
    # Constructor
    ##################################################

    def __init__(self, path):
        '''Initialization of the ImagingReduction instances

        Parameters
        ----------
        path : str
            Path to the directory containing the dataset

        '''

        # expand path
        path = Path(path).expanduser().resolve()

        # zeroth-order reduction validation
        raw = path / 'raw'
        if not raw.exists():
            raise ValueError('No raw/ subdirectory. {0} is not a valid reduction path!'.format(path))

        # init path and name
        self._path = utils.ReductionPath(path)
        self._instrument = 'IRDIS'

        # instrument mode
        self._mode = 'Unknown'

        # configuration
        configfile = Path(vltpf.__file__).parent / 'instruments' / '{}.ini'.format(self._instrument)
        config = configparser.ConfigParser()
        try:
            config.read(configfile)

            # instrument
            self._pixel = float(config.get('instrument', 'pixel'))
            self._nwave = 2

            # calibration
            self._wave_cal_lasers = np.array(eval(config.get('calibration', 'wave_cal_lasers')))
          
            # imaging calibration
            self._default_center = np.array(eval(config.get('calibration-imaging', 'default_center')))
            self._orientation_offset = eval(config.get('calibration-imaging', 'orientation_offset'))

            # reduction parameters
            self._config = {}
            for group in ['reduction', 'reduction-imaging']:
                items = dict(config.items(group))
                self._config.update(items)
                for key, value in items.items():
                    try:
                        val = eval(value)
                    except NameError:
                        val = value
                    self._config[key] = val
        except configparser.Error as e:
            raise ValueError('Error reading configuration file for instrument {0}: {1}'.format(self._instrument, e.message))

        # execution of recipes
        self._recipe_execution = {
            'sort_files': False,
            'sort_frames': False,
            'check_files_association': False,
            'sph_ird_cal_dark': False,
            'sph_ird_cal_detector_flat': False,
            'sph_ird_preprocess_science': False,
            'sph_ird_star_center': False,
            'sph_ird_combine_data': False
        }

        # reload any existing data frames
        self.read_info()

    ##################################################
    # Representation
    ##################################################

    def __repr__(self):
        return '<ImagingReduction, instrument={}, mode={}, path={}>'.format(self._instrument, self._mode, self._path)

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
    def recipe_execution(self):
        return self._recipe_execution

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
        Create static calibrations with esorex
        '''

        config = self._config

        self.sph_ird_cal_dark(silent=config['misc_silent_esorex'])
        self.sph_ird_cal_detector_flat(silent=config['misc_silent_esorex'])


    def preprocess_science(self):
        '''
        Clean and collapse images
        '''

        config = self._config

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

        config = self._config

        self.sph_ird_star_center(high_pass=config['center_high_pass'],
                                 offset=config['center_offset'],
                                 plot=config['misc_plot'])
        self.sph_ird_combine_data(cpix=config['combine_cpix'],
                                  psf_dim=config['combine_psf_dim'],
                                  science_dim=config['combine_science_dim'],
                                  correct_anamorphism=config['combine_correct_anamorphism'],
                                  manual_center=config['combine_manual_center'],
                                  coarse_centering=config['combine_coarse_centering'],
                                  shift_method=config['combine_shift_method'],
                                  save_scaled=config['combine_save_scaled'])


    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''

        config = self._config

        if config['clean']:
            self.sph_ird_clean(delete_raw=config['clean_delete_raw'],
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
    # SPHERE/IRDIS methods
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
        fname = path.preproc / 'files.csv'
        if fname.exists():
            files_info = pd.read_csv(fname, index_col=0)

            # convert times
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

            # update recipe execution
            self._recipe_execution['sort_files'] = True
            if np.any(files_info['PRO CATG'] == 'IRD_MASTER_DARK'):
                self._recipe_execution['sph_ird_cal_dark'] = True
            if np.any(files_info['PRO CATG'] == 'IRD_FLAT_FIELD'):
                self._recipe_execution['sph_ird_cal_detector_flat'] = True

            # update instrument mode
            self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]
        else:
            files_info = None

        fname = path.preproc / 'frames.csv'
        if fname.exists():
            frames_info = pd.read_csv(fname, index_col=(0, 1))

            # convert times
            frames_info['DATE-OBS'] = pd.to_datetime(frames_info['DATE-OBS'], utc=False)
            frames_info['DATE'] = pd.to_datetime(frames_info['DATE'], utc=False)
            frames_info['DET FRAM UTC'] = pd.to_datetime(frames_info['DET FRAM UTC'], utc=False)
            frames_info['TIME START'] = pd.to_datetime(frames_info['TIME START'], utc=False)
            frames_info['TIME'] = pd.to_datetime(frames_info['TIME'], utc=False)
            frames_info['TIME END'] = pd.to_datetime(frames_info['TIME END'], utc=False)

            # update recipe execution
            self._recipe_execution['sort_frames'] = True
        else:
            frames_info = None

        fname = path.preproc / 'frames_preproc.csv'
        if fname.exists():
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
            done = True
            files = frames_info_preproc.index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                file = list(path.preproc.glob('{}.fits'.format(fname)))
                done = done and (len(file) == 1)
            self._recipe_execution['sph_ird_preprocess_science'] = done

            done = True
            files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
                                        (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
            for file, idx in files:
                fname = '{0}_DIT{1:03d}_preproc_centers'.format(file, idx)
                file = list(path.preproc.glob('{}.fits'.format(fname)))
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
        files = path.raw.glob('*.fits')
        files = [f.stem for f in files]

        if len(files) == 0:
            raise ValueError('No raw FITS files in reduction path')

        print(' * found {0} FITS files in {1}'.format(len(files), path.raw))

        # read list of keywords
        keywords = []
        file = open(Path(vltpf.__file__).parent / 'instruments' / 'keywords.dat', 'r')
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
            hdu = fits.open(path.raw / '{}.fits'.format(f))
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
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

        # update instrument mode
        self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'][0]

        # sort by acquisition time
        files_info.sort_values(by='DATE-OBS', inplace=True)

        # save files_info
        files_info.to_csv(path.preproc / 'files.csv')
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

        # raise error when no science frames are present
        if len(sci_files) == 0:
            raise ValueError('This dataset contains no science frame. There should be at least one!')

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
        frames_info.to_csv(path.preproc / 'frames.csv')
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

        posang   = cinfo['INS4 DROT2 POSANG'].unique()

        date = str(cinfo['DATE'][0])[0:10]

        print(' * Object:      {0}'.format(cinfo['OBJECT'][0]))
        print(' * RA / DEC:    {0} / {1}'.format(RA, DEC))
        print(' * Date:        {0}'.format(date))
        print(' * Instrument:  {0}'.format(cinfo['SEQ ARM'][0]))
        print(' * Derotator:   {0}'.format(cinfo['INS4 DROT2 MODE'][0]))
        print(' * Coronagraph: {0}'.format(cinfo['INS COMB ICOR'][0]))
        print(' * Mode:        {0}'.format(cinfo['INS1 MODE'][0]))
        print(' * Filter:      {0}'.format(cinfo['INS COMB IFLT'][0]))
        print(' * DIT:         {0:.2f} sec'.format(cinfo['DET SEQ1 DIT'][0]))
        print(' * NDIT:        {0:.0f}'.format(cinfo['DET NDIT'][0]))
        print(' * Texp:        {0:.2f} min'.format(cinfo['DET SEQ1 DIT'].sum()/60))
        print(' * PA:          {0:.2f}° ==> {1:.2f}° = {2:.2f}°'.format(pa_start, pa_end, np.abs(pa_end-pa_start)))
        print(' * POSANG:      {0}'.format(', '.join(['{:.2f}°'.format(p) for p in posang])))


    def check_files_association(self):
        '''
        Performs the calibration files association as a sanity check.

        Warnings and errors are reported at the end. Execution is
        interupted in case of error.
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'check_files_association', self.recipe_requirements)

        print('Performing file association for calibrations')

        # parameters
        files_info = self._files_info

        # instrument arm
        arm = files_info['SEQ ARM'].unique()
        if len(arm) != 1:
            raise ValueError('Sequence is mixing different instruments: {0}'.format(arm))

        # IRDIS obs mode and filter combination
        modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS1 MODE'].unique()
        if len(modes) != 1:
            raise ValueError('Sequence is mixing different types of observations: {0}'.format(modes))

        filter_combs = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS COMB IFLT'].unique()
        if len(filter_combs) != 1:
            raise ValueError('Sequence is mixing different types of filters combinations: {0}'.format(filter_combs))
        filter_comb = filter_combs[0]

        # specific data frame for calibrations
        # keep static calibrations and sky backgrounds
        calibs = files_info[(files_info['DPR CATG'] == 'CALIB') |
                            ((files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == 'SKY'))]

        ###############################################
        # static calibrations not dependent on DIT
        ###############################################
        error_flag = 0
        warning_flag = 0

        # flat
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS COMB IFLT'] == filter_comb)]
        if len(cfiles) <= 1:
            error_flag += 1
            print(' * Error: there should be more than 1 flat in filter combination {0}'.format(filter_comb))

        ##################################################
        # static calibrations that depend on science DIT
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
                      'It is *highly recommended* to include one to obtain the best data reduction. ' +
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
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_cal_dark', self.recipe_requirements)

        print('Creating darks and backgrounds')

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

                    print(' * {0} in filter {1} with DIT={2:.2f} sec ({3} files)'.format(ctype, cfilt, DIT, len(cfiles)))

                    # create sof
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
                            '--ird.master_dark.outfilename={0}/{1}.fits'.format(path.calib, dark_file),
                            '--ird.master_dark.badpixfilename={0}/{1}.fits'.format(path.calib, bpm_file),
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
        files_info.to_csv(path.preproc / 'files.csv')

        # update recipe execution
        self._recipe_execution['sph_ird_cal_dark'] = True


    def sph_ird_cal_detector_flat(self, silent=True):
        '''
        Create the detector flat calibrations

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_cal_detector_flat', self.recipe_requirements)

        print('Creating flats')

        # parameters
        path = self._path
        files_info = self._files_info

        # get list of files
        calibs = files_info[np.logical_not(files_info['PROCESSED']) &
                            ((files_info['DPR TYPE'] == 'FLAT,LAMP') |
                             (files_info['DPR TECH'] == 'IMAGE'))]
        filter_combs = calibs['INS COMB IFLT'].unique()

        for cfilt in filter_combs:
            cfiles = calibs[calibs['INS COMB IFLT'] == cfilt]
            files = cfiles.index

            print(' * filter {0} ({1} files)'.format(cfilt, len(cfiles)))

            # create sof
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
        files_info.to_csv(path.preproc / 'files.csv')

        # update recipe execution
        self._recipe_execution['sph_ird_cal_detector_flat'] = True


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

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_preprocess_science', self.recipe_requirements)

        print('Pre-processing science files')

        # parameters
        path = self._path
        files_info = self._files_info
        frames_info = self._frames_info

        # clean before we start
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

            bpm = toolbox.compute_bad_pixel_map(bpm_files)

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
            raise ValueError('There should be exactly 1 flat file. Found {0}.'.format(len(flat_file)))
        flat = fits.getdata(path.calib / '{}.fits'.format(flat_file.index[0]))

        # final dataframe
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

                    if len(dfiles) == 0:
                        # issue a warning if absolutely no background is found
                        print('Warning: no background has been found. Pre-processing will continue but data quality will likely be affected')
                        bkg = np.zeros((1024, 2048))
                    elif len(dfiles) == 1:
                        bkg = fits.getdata(path.calib / '{}.fits'.format(dfiles.index[0]))
                    elif len(dfiles) > 1:
                        # FIXME: handle cases when multiple backgrounds are found?
                        raise ValueError('Unexpected number of background files ({0})'.format(len(dfiles)))

                # process files
                for idx, (fname, finfo) in enumerate(sfiles.iterrows()):
                    # frames_info extract
                    finfo = frames_info.loc[(fname, slice(None)), :]

                    print(' * file {0}/{1}: {2}, NDIT={3}'.format(idx+1, len(sfiles), fname, len(finfo)))

                    # read data
                    print('   ==> read data')
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
                                nimg = np.empty((NDIT_new, 1024, 2048), dtype=img.dtype)
                                for f in range(NDIT_new):
                                    nimg[f] = np.mean(img[f*coadd_value:(f+1)*coadd_value], axis=0)
                                img = nimg

                                frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'coadd', coadd_value=coadd_value)
                            else:
                                raise ValueError('Unknown collapse type {0}'.format(collapse_type))
                        else:
                            frames_info_new = toolbox.collapse_frames_info(finfo, fname, 'none')

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
                            frame = imutils.fix_badpix(frame, bpm, npix=12, weight=True)

                            # additional sigma clipping to remove transitory bad pixels
                            # not done for OBJECT,FLUX because PSF peak can be clipped
                            if (typ != 'OBJECT,FLUX'):
                                frame = imutils.sigma_filter(frame, box=7, nsigma=4, iterate=False)

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
                        fits.writeto(path.preproc / '{}_DIT{:03d}_preproc.fits'.format(fname, f), frame, hdr,
                                     overwrite=True, output_verify='silentfix')

                    print()

            print()

        # sort and save final dataframe
        frames_info_preproc.sort_values(by='TIME', inplace=True)
        frames_info_preproc.to_csv(path.preproc / 'frames_preproc.csv')

        self._frames_info_preproc = frames_info_preproc

        # update recipe execution
        self._recipe_execution['sph_ird_preprocess_science'] = True


    def sph_ird_star_center(self, high_pass=False, offset=(0, 0), plot=True):
        '''Determines the star center for all frames where a center can be
        determined (OBJECT,CENTER and OBJECT,FLUX)

        Parameters
        ----------
        high_pass : bool
            Apply high-pass filter to the image before searching for the satelitte spots.
            Default is False

        offset : tuple
            Apply an (x,y) offset to the default center position, for the waffle centering.
            The offset will move the search box of the waffle spots by the amount of
            specified pixels in each direction. Default is no offset

        plot : bool
            Display and save diagnostic plot for quality check. Default is True

        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_star_center', self.recipe_requirements)

        print('Star centers determination')

        # parameters
        path = self._path
        pixel = self._pixel
        orientation_offset = self._orientation_offset
        center_guess = self._default_center
        frames_info = self._frames_info_preproc        

        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)
        wave = np.array(wave)

        # start with OBJECT,FLUX
        flux_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,FLUX']
        if len(flux_files) != 0:
            for file, idx in flux_files.index:
                print('  ==> OBJECT,FLUX: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube, hdr = fits.getdata(path.preproc / '{}.fits'.format(fname), header=True)

                # centers
                if plot:
                    save_path = path.products / '{}_PSF_fitting.pdf'.format(fname)
                else:
                    save_path = None
                img_center = toolbox.star_centers_from_PSF_img_cube(cube, wave, pixel, save_path=save_path)

                # save
                fits.writeto(path.preproc / '{}_centers.fits'.format(fname), img_center, overwrite=True)
                print()

        # then OBJECT,CENTER
        starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
        if len(starcen_files) != 0:
            for file, idx in starcen_files.index:
                print('  ==> OBJECT,CENTER: {0}'.format(file))

                # read data
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube, hdr = fits.getdata(path.preproc / '{}.fits'.format(fname), header=True)

                # coronagraph
                coro_name = starcen_files.loc[(file, idx), 'INS COMB ICOR']
                if coro_name == 'N_NS_CLEAR':
                    coro = False
                else:
                    coro = True

                # centers
                waffle_orientation = hdr['HIERARCH ESO OCS WAFFLE ORIENT']
                if plot:
                    save_path = path.products / '{}_spots_fitting.pdf'.format(fname)
                else:
                    save_path = None
                spot_center, spot_dist, img_center \
                    = toolbox.star_centers_from_waffle_img_cube(cube, wave, waffle_orientation, center_guess,
                                                                pixel, orientation_offset, high_pass=high_pass, 
                                                                center_offset=offset, coro=coro, save_path=save_path)

                # save
                fits.writeto(path.preproc / '{}_centers.fits'.format(fname), img_center, overwrite=True)
                print()

        # update recipe execution
        self._recipe_execution['sph_ird_star_center'] = True


    def sph_ird_combine_data(self, cpix=True, psf_dim=80, science_dim=290, correct_anamorphism=True,
                             shift_method='fft', manual_center=None, coarse_centering=False, save_scaled=False):
        '''Combine and save the science data into final cubes

        All types of data are combined independently: PSFs
        (OBJECT,FLUX), star centers (OBJECT,CENTER) and standard
        coronagraphic images (OBJECT). For each type of data, the
        method saves 4 or 5 different files:

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
        pre-defined center that a representative of the typical center
        of the coronagraph.

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
            frames. This should be an array of either 2 or nwavex2
            values. For OBJECT,FLUX frames, the PSF is always
            recentered. Default is None

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
            to True when coarse_centering is set to True. The default
            is False
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_combine_data', self.recipe_requirements)

        print('Combine science data')

        # parameters
        path = self._path
        nwave = self._nwave
        frames_info = self._frames_info_preproc

        # wavelength
        filter_comb = frames_info['INS COMB IFLT'].unique()[0]
        wave, bandwidth = transmission.wavelength_bandwidth_filter(filter_comb)
        wave = np.array(wave)

        fits.writeto(path.products / 'wavelength.fits', wave, overwrite=True)

        # max images size
        if psf_dim > 1024:
            print('Warning: psf_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            psf_dim = 1024

        if science_dim > 1024:
            print('Warning: science_dim cannot be larger than 1024 pix. A value of 1024 will be used.')
            science_dim = 1024

        # centering configuration
        if coarse_centering:
            print('Warning: images will be coarsely centered without any interpolation. Automatic settings for coarse centering: shift_method=\'roll\', cpix=True, correct_anamorphism=False, save_scaled=False')
            shift_method = 'roll'
            cpix = True
            correct_anamorphism = False
            save_scaled = False

        if manual_center is not None:
            manual_center = np.array(manual_center)
            
            if (manual_center.shape != (2,)) and (manual_center.shape != (nwave, 2)):
                raise ValueError('manual_center does not have the right number of dimensions.')

            if manual_center.shape == (2,):
                manual_center = np.full((nwave, 2), manual_center)

            print('Warning: images will be centered using the user-provided center ({},{})'.format(*manual_center[0]))

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
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube = fits.getdata(path.preproc / '{}.fits'.format(fname))
                
                cfile = path.preproc / '{}_centers.fits'.format(fname)
                if cfile.exists():
                    centers = fits.getdata(cfile)
                else:
                    print('Warning: sph_ifs_star_center() has not been executed. Images will be centered using default center ({},{})'.format(*self._default_center))
                    centers = np.full((nwave, 2), self._default_center)

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)

                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                psf_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                psf_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    psf_cube[wave_idx, file_idx] = nimg[:psf_dim, :psf_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = psf_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        psf_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = psf_cube[wave_idx, file_idx]
                        psf_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            flux_files.to_csv(path.products / 'psf_frames.csv')
            fits.writeto(path.products / 'psf_cube.fits', psf_cube, overwrite=True)
            fits.writeto(path.products / 'psf_parang.fits', psf_parang, overwrite=True)
            fits.writeto(path.products / 'psf_derot.fits', psf_derot, overwrite=True)
            if save_scaled:
                fits.writeto(path.products / 'psf_cube_scaled.fits', psf_cube_scaled, overwrite=True)

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
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                cube = fits.getdata(path.preproc / '{}.fits'.format(fname))

                # use manual center if explicitely requested
                if manual_center is not None:
                    centers = manual_center
                else:
                    # otherwise read center data
                    centers = fits.getdata(path.preproc / '{}_centers.fits'.format(fname))
                
                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    centers = centers.astype(np.int)
                
                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                cen_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                cen_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    cen_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = cen_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        cen_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = cen_cube[wave_idx, file_idx]
                        cen_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            starcen_files.to_csv(path.products / 'starcenter_frames.csv')
            fits.writeto(path.products / 'starcenter_cube.fits', cen_cube, overwrite=True)
            fits.writeto(path.products / 'starcenter_parang.fits', cen_parang, overwrite=True)
            fits.writeto(path.products / 'starcenter_derot.fits', cen_derot, overwrite=True)
            if save_scaled:
                fits.writeto(path.products / 'starcenter_cube_scaled.fits', cen_cube_scaled, overwrite=True)

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

            # use manual center if explicitely requested
            if manual_center is not None:
                centers = manual_center
            else:
                # otherwise, look whether we have an OBJECT,CENTER frame
            
                # FIXME: ticket #12. Use first DIT of first OBJECT,CENTER
                # in the sequence, but it would be better to be able to
                # select which CENTER to use
                starcen_files = frames_info[frames_info['DPR TYPE'] == 'OBJECT,CENTER']
                if len(starcen_files) == 0:
                    print('Warning: no OBJECT,CENTER file in the dataset. Images will be centered using default center ({},{})'.format(*self._default_center))
                    centers = np.full((nwave, 2), self._default_center)
                    
                    # null value for Dithering Motion Stage
                    dms_dx_ref = 0
                    dms_dy_ref = 0
                else:
                    fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                    centers = fits.getdata(path.preproc / fname)
                    fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                    
                    # Dithering Motion Stage for star center: value is in micron,
                    # and the pixel size is 18 micron
                    dms_dx_ref = starcen_files['INS1 PAC X'][0] / 18
                    dms_dy_ref = starcen_files['INS1 PAC Y'][0] / 18

            # make sure we have only integers if user wants coarse centering
            if coarse_centering:
                centers = centers.astype(np.int)

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
                fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
                files = list(path.preproc.glob('{}*.fits'.format(fname)))
                cube = fits.getdata(files[0])

                # neutral density
                ND = frames_info.loc[(file, idx), 'INS4 FILT2 NAME']
                w, attenuation = transmission.transmission_nd(ND, wave=wave)

                # DIT, angles, etc
                DIT = frames_info.loc[(file, idx), 'DET SEQ1 DIT']
                sci_parang[file_idx] = frames_info.loc[(file, idx), 'PARANG']
                sci_derot[file_idx] = frames_info.loc[(file, idx), 'DEROT ANGLE']

                # Dithering Motion Stage for star center: value is in micron,
                # and the pixel size is 18 micron
                dms_dx = frames_info.loc[(file, idx), 'INS1 PAC X'] / 18
                dms_dy = frames_info.loc[(file, idx), 'INS1 PAC Y'] / 18

                # make sure we have only integers if user wants coarse centering
                if coarse_centering:
                    dms_dx = np.int(dms_dx)
                    dms_dy = np.int(dms_dy)
                
                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    # DMS contribution
                    cx = cx + dms_dx_ref + dms_dx
                    cy = cy + dms_dy_ref + dms_dy

                    print(cx, cy)
                    
                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method=shift_method)
                    nimg = nimg / DIT / attenuation[wave_idx]

                    sci_cube[wave_idx, file_idx] = nimg[:science_dim, :science_dim]

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = sci_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        sci_cube[wave_idx, file_idx] = nimg

                    # wavelength-scaled version
                    if save_scaled:
                        nimg = sci_cube[wave_idx, file_idx]
                        sci_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method=shift_method)

            # save final cubes
            object_files.to_csv(path.products / 'science_frames.csv')
            fits.writeto(path.products / 'science_cube.fits', sci_cube, overwrite=True)
            fits.writeto(path.products / 'science_parang.fits', sci_parang, overwrite=True)
            fits.writeto(path.products / 'science_derot.fits', sci_derot, overwrite=True)
            if save_scaled:
                fits.writeto(path.products / 'science_cube_scaled.fits', sci_cube_scaled, overwrite=True)

            # delete big cubes
            del sci_cube
            if save_scaled:
                del sci_cube_scaled

            print()

        # update recipe execution
        self._recipe_execution['sph_ird_combine_data'] = True


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
        if path.tmp.exists():
            shutil.rmtree(path.tmp, ignore_errors=True)

        # sof
        if path.sof.exists():
            shutil.rmtree(path.sof, ignore_errors=True)

        # calib
        if path.calib.exists():
            shutil.rmtree(path.calib, ignore_errors=True)

        # preproc
        if path.preproc.exists():
            shutil.rmtree(path.preproc, ignore_errors=True)

        # raw
        if delete_raw:
            if path.raw.exists():
                shutil.rmtree(path.raw, ignore_errors=True)

        # products
        if delete_products:
            if path.products.exists():
                shutil.rmtree(path.products, ignore_errors=True)
