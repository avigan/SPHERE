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


class ImagingReduction(object):
    '''SPHERE/IRDIS imaging reduction object. It handles the dual-band
    imaging (DBI) and classifcal imaging (CI) observing modes.

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
        'sph_ird_preprocess_science': ['sort_files', 'sort_frames', 'sph_ird_cal_dark', 'sph_ird_cal_detector_flat'],
        'sph_ird_star_center': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science'],
        'sph_ird_combine_data': ['sort_files', 'sort_frames', 'sph_ird_preprocess_science', 'sph_ird_star_center']
    }
    
    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, path):
        '''Initialization of the ImagingReduction instances

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
        self._instrument = 'IRDIS'
        
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
            'check_files_association': False
        }
        
        # reload any existing data frames
        self.read_info()
    
    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return '<ImagingReduction, instrument={0}, path={1}>'.format(self._instrument, self._path)
    
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

        # sort files and frames
        self.sort_files()
        self.sort_frames()

        # sanity check
        self.check_files_association()
        
    
    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        config = self._config
        
        self.sph_ird_cal_dark(silent=config['silent'])
        self.sph_ird_cal_detector_flat(silent=config['silent'])

    
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
                                 display=config['center_display'],
                                 save=config['center_save'])
        self.sph_ird_combine_data(cpix=config['combine_cpix'],
                                  psf_dim=config['combine_psf_dim'],
                                  science_dim=config['combine_science_dim'],
                                  correct_anamorphism=config['combine_correct_anamorphism'],
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
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_cal_dark', self.recipe_requirements)
        
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
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

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
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

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
        '''Determines the star center for all frames where a center can be
        determined (OBJECT,CENTER and OBJECT,FLUX)

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
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_star_center', self.recipe_requirements)
        
        print('Star centers determination')

        # parameters
        path = self._path
        pixel = self._pixel
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
                img_center = toolbox.star_centers_from_PSF_cube(cube, wave, pixel, display=display, save_path=save_path)

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
                spot_center, spot_dist, img_center \
                    = toolbox.star_centers_from_waffle_cube(cube, wave, 'IRDIS', waffle_orientation,
                                                            high_pass=high_pass, display=display,
                                                            save_path=save_path)

                # save
                fits.writeto(os.path.join(path.preproc, fname+'_centers.fits'), img_center, overwrite=True)
                print()

        # update recipe execution
        self._recipe_execution['sph_ird_star_center'] = True


    def sph_ird_combine_data(self, cpix=True, psf_dim=80, science_dim=290, correct_anamorphism=True, save_scaled=False):
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

        save_scaled : bool    
            Also save the wavelength-rescaled cubes. Makes the process
            much longer. The default is False

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

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = psf_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        psf_cube[wave_idx, file_idx] = nimg
                    
                    # wavelength-scaled version
                    if save_scaled:
                        nimg = psf_cube[wave_idx, file_idx]
                        psf_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

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

                    # correct anamorphism
                    if correct_anamorphism:
                        nimg = cen_cube[wave_idx, file_idx]
                        nimg = imutils.scale(nimg, (1.0000, 1.0062), method='interp')
                        cen_cube[wave_idx, file_idx] = nimg
                    
                    # wavelength-scaled version
                    if save_scaled:
                        nimg = cen_cube[wave_idx, file_idx]
                        cen_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

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
                print(' ==> no OBJECT,CENTER file in the data set. Images cannot be accurately centred. ' +
                      'They will just be combined.')

                centers = np.full((nwave, 2), cc)

                # null value for Dithering Motion Stage
                dms_dx_ref = 0
                dms_dy_ref = 0
            else:
                fname = '{0}_DIT{1:03d}_preproc_centers.fits'.format(starcen_files.index.values[0][0], starcen_files.index.values[0][1])
                centers = fits.getdata(os.path.join(path.preproc, fname))

                # Dithering Motion Stage for star center: value is in micron,
                # and the pixel size is 18 micron
                dms_dx_ref = starcen_files['INS1 PAC X'][0] / 18
                dms_dy_ref = starcen_files['INS1 PAC Y'][0] / 18
                
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
                print('  ==> file {0}/{1}: {2}, DIT={3}'.format(file_idx+1, len(object_files), file, idx))

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

                # Dithering Motion Stage for star center: value is in micron,
                # and the pixel size is 18 micron
                dms_dx = frames_info.loc[(file, idx), 'INS1 PAC X'] / 18
                dms_dy = frames_info.loc[(file, idx), 'INS1 PAC Y'] / 18
                
                # center frames
                for wave_idx, img in enumerate(cube):
                    cx, cy = centers[wave_idx, :]

                    # DMS contribution
                    cx = cx + dms_dx_ref + dms_dx
                    cy = cy + dms_dy_ref + dms_dy

                    img  = img.astype(np.float)
                    nimg = imutils.shift(img, (cc-cx, cc-cy), method='fft')
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
                        sci_cube_scaled[wave_idx, file_idx] = imutils.scale(nimg, wave[0]/wave[wave_idx], method='fft')

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
