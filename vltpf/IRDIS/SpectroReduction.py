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

import vltpf
import vltpf.utils.imutils as imutils
import vltpf.utils.aperture as aperture
import vltpf.transmission as transmission
import vltpf.ReductionPath as ReductionPath
import vltpf.toolbox as toolbox


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
        'sph_ird_wave_calib': ['sort_files']
    }
    
    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, path):
        '''Initialization of the SpectroReduction instances

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
        package_directory = os.path.dirname(os.path.abspath(vltpf.__file__))
        configfile = os.path.join(package_directory, 'instruments', self._instrument+'.ini')
        config = configparser.ConfigParser()
        try:
            config.read(configfile)

            # instrument
            self._pixel = float(config.get('instrument', 'pixel'))
            self._nwave = int(config.get('instrument', 'nwave'))
            self._wave_cal_lasers = [float(w) for w in config.get('calibration', 'wave_cal_lasers').split(',')]

            # reduction
            self._config = dict(config.items('reduction-spectro'))
            for key, value in self._config.items():
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
            'sph_ifs_cal_dark': False,
            'sph_ifs_cal_detector_flat': False,
            'sph_ird_wave_calib': False
        }
        
        # reload any existing data frames
        self.read_info()
    
    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return '<SpectroReduction, instrument={0}, path={1}>'.format(self._instrument, self._path)
    
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
        Create static calibrations with esorex
        '''

        config = self._config
        
        self.sph_ird_cal_dark(silent=config['silent'])
        self.sph_ird_cal_detector_flat(silent=config['silent'])
        self.sph_ird_wave_calib(silent=config['silent'])

    
    def preprocess_science(self):
        '''
        Clean and collapse images
        '''
        
        config = self._config
        
        pass
    

    def process_science(self):
        '''
        Perform star center, combine cubes into final (x,y,time,lambda)
        cubes, correct anamorphism and scale the images
        '''
        
        config = self._config
        
        pass
    
    
    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''
        
        config = self._config

        pass
    
        
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
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)
            
            # update recipe execution
            self._recipe_execution['sort_files'] = True
            # if np.any(files_info['PRO CATG'] == 'IRD_MASTER_DARK'):
            #     self._recipe_execution['sph_ird_cal_dark'] = True
            # if np.any(files_info['PRO CATG'] == 'IRD_FLAT_FIELD'):
            #     self._recipe_execution['sph_ird_cal_detector_flat'] = True
        else:
            files_info = None

        fname = os.path.join(path.preproc, 'frames.csv')
        if os.path.exists(fname):
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

        fname = os.path.join(path.preproc, 'frames_preproc.csv')
        if os.path.exists(fname):
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
        # if frames_info_preproc is not None:
        #     done = True
        #     files = frames_info_preproc.index
        #     for file, idx in files:
        #         fname = '{0}_DIT{1:03d}_preproc'.format(file, idx)
        #         file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
        #         done = done and (len(file) == 1)
        #     self._recipe_execution['sph_ird_preprocess_science'] = done

        #     done = True
        #     files = frames_info_preproc[(frames_info_preproc['DPR TYPE'] == 'OBJECT,FLUX') |
        #                                 (frames_info_preproc['DPR TYPE'] == 'OBJECT,CENTER')].index
        #     for file, idx in files:
        #         fname = '{0}_DIT{1:03d}_preproc_centers'.format(file, idx)
        #         file = glob.glob(os.path.join(path.preproc, fname+'.fits'))
        #         done = done and (len(file) == 1)
        #     self._recipe_execution['sph_ird_star_center'] = done

        
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
        package_directory = os.path.dirname(os.path.abspath(vltpf.__file__))
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
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=False)

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
        path = self._path
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
        if filter_comb == 'S_LR':
            mode_short = 'LRS'
        elif filter_comb == 'S_MR':
            mode_short = 'MRS'
        else:
            raise ValueError('Unknown IRDIS-LSS mode {0}'.format(filter_comb))

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
        
        # wave
        cfiles = calibs[(calibs['DPR TYPE'] == 'LAMP,WAVE') & (calibs['INS COMB IFLT'] == filter_comb)]
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
            files_info.drop(time_delta[1:].index, inplace=True)
        
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
                      'usually provide a cleaner data reduction')

        # error reporting
        print('There are {0} warning(s) and {1} error(s) in the classification of files'.format(warning_flag, error_flag))
        if error_flag:
            raise ValueError('There is {0} errors that should be solved before proceeding'.format(error_flag))

        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))
        self._files_info = files_info

    
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
        calibs = files_info[np.logical_not(files_info['PROCESSED']) &
                            (files_info['DPR TYPE'] == 'FLAT,LAMP')]
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

    
    def sph_ird_wave_calib(self, silent=True):
        '''
        Create the wavelength calibration

        Parameters
        ----------
        silent : bool
            Suppress esorex output. Default is True
        '''

        # check if recipe can be executed
        toolbox.check_recipe_execution(self._recipe_execution, 'sph_ird_wave_calib', self.recipe_requirements)
        
        print('Creating wavelength calibration')

        # parameters
        path = self._path
        files_info = self._files_info
        
        # get list of files
        wave_file = files_info[np.logical_not(files_info['PROCESSED']) & (files_info['DPR TYPE'] == 'LAMP,WAVE')]
        if len(wave_file) != 1:
            raise ValueError('There should be exactly 1 raw wavelength calibration file. Found {0}.'.format(len(wave_file)))
        
        DIT = wave_file['DET SEQ1 DIT'][0]
        dark_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_MASTER_DARK') & 
                               (files_info['DPR CATG'] == 'CALIB') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
        if len(dark_file) == 0:
            raise ValueError('There should at least 1 dark file for wavelength calibration. Found none.')

        mode = wave_file['INS COMB IFLT'][0]
        flat_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_FLAT_FIELD')]
        if len(flat_file) == 0:
            raise ValueError('There should at least 1 flat file for wavelength calibration. Found none.')
        
        bpm_file = files_info[files_info['PROCESSED'] & (files_info['PRO CATG'] == 'IRD_NON_LINEAR_BADPIXELMAP')]
        if len(flat_file) == 0:
            raise ValueError('There should at least 1 bad pixel map file for wavelength calibration. Found none.')
        
        # products
        wav_file = 'wave_calib'
        
        # esorex parameters
        if mode == 'S_LR':
            # create standard sof in LRS
            sof = os.path.join(path.sof, 'wave.sof')
            file = open(sof, 'w')
            file.write('{0}{1}.fits     {2}\n'.format(path.raw, wave_file, 'IRD_WAVECALIB_RAW'))
            file.write('{0}{1}.fits     {2}\n'.format(path.calib, dark_file.index[0], 'IRD_MASTER_DARK'))
            file.write('{0}{1}.fits     {2}\n'.format(path.calib, flat_file.index[0], 'IRD_FLAT_FIELD'))
            file.write('{0}{1}.fits     {2}\n'.format(path.calib, bpm_file.index[0], 'IRD_STATIC_BADPIXELMAP'))
            file.close()
            
            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_wave_calib',
                    '--ird.wave_calib.column_width=200',
                    '--ird.wave_calib.grism_mode=FALSE',
                    '--ird.wave_calib.threshold=2000',
                    '--ird.wave_calib.number_lines=6',
                    '--ird.wave_calib.outfilename={0}{1}.fits'.format(path.calib, wav_file),
                    sof]
        elif mode == 'S_MR':            
            # masking of second order spectrum in MRS
            wave_fname = wave_file.index[0]
            wave_data, hdr = fits.getdata(os.path.join(path.raw, wave_fname+'.fits'), header=True)
            wave_data = wave_data.squeeze()
            wave_data[:60, :] = 0
            fits.writeto(os.path.join(path.preproc, wave_fname+'_masked.fits'), wave_data, hdr, overwrite=True, 
                         output_verify='silentfix')
            
            # create sof using the masked file
            sof = os.path.join(path.sof, 'wave.sof')
            file = open(sof, 'w')
            file.write('{0}{1}_masked.fits {2}\n'.format(path.preproc, wave_fname, 'IRD_WAVECALIB_RAW'))
            file.write('{0}{1}.fits        {2}\n'.format(path.calib, dark_file.index[0], 'IRD_MASTER_DARK'))
            file.write('{0}{1}.fits        {2}\n'.format(path.calib, flat_file.index[0], 'IRD_FLAT_FIELD'))
            file.write('{0}{1}.fits        {2}\n'.format(path.calib, bpm_file.index[0], 'IRD_STATIC_BADPIXELMAP'))
            file.close()

            args = ['esorex',
                    '--no-checksum=TRUE',
                    '--no-datamd5=TRUE',
                    'sph_ird_wave_calib',
                    '--ird.wave_calib.column_width=200',
                    '--ird.wave_calib.grism_mode=TRUE',
                    '--ird.wave_calib.threshold=1000',
                    '--ird.wave_calib.number_lines=5',
                    '--ird.wave_calib.outfilename={0}{1}.fits'.format(path.calib, wav_file),
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
        files_info.loc[wav_file, 'INS COMB IFLT'] = wave_file['INS COMB IFLT'][0]
        files_info.loc[wav_file, 'INS4 FILT2 NAME'] = wave_file['INS4 FILT2 NAME'][0]
        files_info.loc[wav_file, 'INS1 MODE'] = wave_file['INS1 MODE'][0]
        files_info.loc[wav_file, 'INS1 FILT NAME'] = wave_file['INS1 FILT NAME'][0]
        files_info.loc[wav_file, 'INS1 OPTI2 NAME'] = wave_file['INS1 OPTI2 NAME'][0]
        files_info.loc[wav_file, 'DET SEQ1 DIT'] = wave_file['DET SEQ1 DIT'][0]
        files_info.loc[wav_file, 'PROCESSED'] = True
        files_info.loc[wav_file, 'PRO CATG'] = 'IRD_WAVECALIB'
        
        # save
        files_info.to_csv(os.path.join(path.preproc, 'files.csv'))

        # update recipe execution
        self._recipe_execution['sph_ird_wave_calib'] = True
