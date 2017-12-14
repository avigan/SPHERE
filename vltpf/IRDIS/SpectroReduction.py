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
        'check_files_association': ['sort_files']
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
                self._config[key] = eval(value)
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

        # make sure we have sub-directories
        self._path.create_subdirectories()
                
        pass
    
    
    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        config = self._config
        
        pass

    
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
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
            files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)
            
            # update recipe execution
            # self._recipe_execution['sort_files'] = True
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
