import pandas as pd
import logging
import numpy as np
import collections
import configparser
import shutil
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import requests
import io

from astropy.io import fits
from astropy.time import Time
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages

import sphere
import sphere.utils as utils
import sphere.toolbox as toolbox

_log = logging.getLogger(__name__)

# WFS wavelength
wave_wfs = 500e-9


class Reduction(object):
    '''
    SPHERE/SPARTA dataset reduction class

    The analysis and plotting code of this class was originally
    developed by Julien Milli (ESO/IPAG) and based on SAXO tools
    from Jean-François Sauvage (ONERA). See:

    https://github.com/jmilou/sparta

    for the code from Julien Milli.
    '''

    ##################################################
    # Class variables
    ##################################################

    # specify for each recipe which other recipes need to have been executed before
    recipe_requirements = collections.OrderedDict([
        ('sort_files', []),
        ('sph_sparta_dtts', ['sort_files']),
        ('sph_sparta_wfs_parameters', ['sort_files']),
        ('sph_sparta_atmospheric_parameters', ['sort_files']),
        ('sph_query_databases', ['sort_files']),
        ('sph_sparta_plot', ['sort_files', 'sph_sparta_dtts', 'sph_sparta_wfs_parameters', 'sph_sparta_atmospheric_parameters']),
        ('sph_ifs_clean', [])
    ])

    ##################################################
    # Constructor
    ##################################################

    def __new__(cls, path, log_level='info', sphere_handler=None):
        '''
        Custom instantiation for the class

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
            reduction = super(Reduction, cls).__new__(cls)

        #
        # basic init
        #

        # init path
        reduction._path = utils.ReductionPath(path)
        
        # instrument and mode
        reduction._instrument = 'SPARTA'

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
        
        reduction._logger.info('Creating SPARTA reduction at path {}'.format(path))

        #
        # configuration
        #
        reduction._logger.debug('> read default configuration')
        configfile = f'{Path(sphere.__file__).parent}/instruments/{reduction._instrument}.ini'
        config = configparser.ConfigParser()

        reduction._logger.debug('Read configuration')
        config.read(configfile)

        # reduction parameters
        reduction._config = dict(config.items('reduction'))
        for key, value in reduction._config.items():
            try:
                val = eval(value)
            except NameError:
                val = value
            reduction._config[key] = val

        #
        # reduction and recipe status
        #
        reduction._status = sphere.INIT
        reduction._recipes_status = collections.OrderedDict()

        for recipe in reduction.recipe_requirements.keys():
            reduction._update_recipe_status(recipe, sphere.NOTSET)
        
        # reload any existing data frames
        reduction._read_info()

        reduction._logger.warning('#########################################################')
        reduction._logger.warning('#                        WARNING!                       #')
        reduction._logger.warning('# Support for SPARTA files is preliminary. The current  #')
        reduction._logger.warning('# format of product files may change in future versions #')
        reduction._logger.warning('# of the pipeline until an appropriate format is found. #')
        reduction._logger.warning('# Please do not blindly rely on the current format.     #')
        reduction._logger.warning('#########################################################')
        
        #
        # return instance
        #
        return reduction

    ##################################################
    # Representation
    ##################################################

    def __repr__(self):
        return '<Reduction, instrument={}, path={}, log={}>'.format(self._instrument, self._path, self.loglevel)

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
    def path(self):
        return self._path

    @property
    def files_info(self):
        return self._files_info

    @property
    def dtts_info(self):
        return self._dtts_info
        
    @property
    def visloop_info(self):
        return self._visloop_info
        
    @property
    def irloop_info(self):
        return self._irloop_info
        
    @property
    def atmospheric_info(self):
        return self._atmos_info
    
    @property
    def recipe_status(self):
        return self._recipes_status

    @property
    def config(self):
        return self._config

    @property
    def status(self):
        return self._status
        
    ##################################################
    # Private methods
    ##################################################

    def _read_info(self):
        '''
        Read the files, calibs and frames information from disk

        files_info : dataframe
            The data frame with all the information on files

        This function is not supposed to be called directly by the user.
        '''

        self._logger.info('Read existing reduction information')
        
        # path
        path = self.path

        # files info
        fname = path.preproc / 'files.csv'
        if fname.exists():
            self._logger.debug('> read files.csv')
            
            files_info = pd.read_csv(fname, index_col=0)

            # convert times
            files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
            files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)

            # update recipe execution
            self._update_recipe_status('sort_files', sphere.SUCCESS)
        else:
            files_info = None

        # DTTS info
        fname = path.products / 'dtts_frames.csv'
        if fname.exists():
            self._logger.debug('> read dtts_frames.csv')
            
            dtts_info = pd.read_csv(fname, index_col=0)

            # convert times
            dtts_info['DATE-OBS'] = pd.to_datetime(dtts_info['DATE-OBS'], utc=False)
            dtts_info['DATE'] = pd.to_datetime(dtts_info['DATE'], utc=False)
            dtts_info['TIME'] = pd.to_datetime(dtts_info['TIME'], utc=False)

            # update recipe execution
            self._update_recipe_status('sph_sparta_dtts', sphere.SUCCESS)
        else:
            dtts_info = None

        # VisLoop info
        fname = path.products / 'visloop_info.csv'
        visloop = False
        if fname.exists():
            self._logger.debug('> read visloop_info.csv')
            
            visloop_info = pd.read_csv(fname, index_col=0)

            # convert times
            visloop_info['DATE-OBS'] = pd.to_datetime(visloop_info['DATE-OBS'], utc=False)
            visloop_info['DATE'] = pd.to_datetime(visloop_info['DATE'], utc=False)
            visloop_info['TIME'] = pd.to_datetime(visloop_info['TIME'], utc=False)

            visloop = True
        else:
            visloop_info = None

        # IRLoop info
        fname = path.products / 'irloop_info.csv'
        irloop = False
        if fname.exists():
            self._logger.debug('> read irloop_info.csv')
            
            irloop_info = pd.read_csv(fname, index_col=0)

            # convert times
            irloop_info['DATE-OBS'] = pd.to_datetime(irloop_info['DATE-OBS'], utc=False)
            irloop_info['DATE'] = pd.to_datetime(irloop_info['DATE'], utc=False)
            irloop_info['TIME'] = pd.to_datetime(irloop_info['TIME'], utc=False)

            irloop = True
        else:
            irloop_info = None

        # update recipe execution
        if visloop and irloop:
            self._update_recipe_status('sph_sparta_wfs_parameters', sphere.SUCCESS)
        else:
            self._update_recipe_status('sph_sparta_wfs_parameters', sphere.NOTSET)

        # Atmospheric info
        fname = path.products / 'atmospheric_info.csv'
        if fname.exists():
            self._logger.debug('> read atmospheric_info.csv')
            
            atmos_info = pd.read_csv(fname, index_col=0)

            # convert times
            atmos_info['DATE-OBS'] = pd.to_datetime(atmos_info['DATE-OBS'], utc=False)
            atmos_info['DATE'] = pd.to_datetime(atmos_info['DATE'], utc=False)
            atmos_info['TIME'] = pd.to_datetime(atmos_info['TIME'], utc=False)

            # update recipe execution
            self._update_recipe_status('sph_sparta_atmospheric_parameters', sphere.SUCCESS)
        else:
            atmos_info = None

        # save data frames in instance variables
        self._files_info   = files_info
        self._dtts_info    = dtts_info
        self._visloop_info = visloop_info
        self._irloop_info  = irloop_info
        self._atmos_info   = atmos_info

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

        self.sort_files()

        
    def create_static_calibrations(self):
        '''
        Create static calibrations
        '''
        
        self._logger.info('====> Static calibrations <====')                
        self._logger.warning('No static calibrations for SPARTA data')


    def preprocess_science(self):
        '''
        Pre-processing of data
        '''

        self._logger.info('====> Science pre-processing <====')
        self._logger.warning('No pre-processing required for SPARTA data')


    def process_science(self):
        '''
        Process the SPARTA files
        '''
        
        self._logger.info('====> Science processing <====')

        config = self._config

        self.sph_sparta_dtts(plot=config['misc_plot'])
        self.sph_sparta_wfs_parameters()
        self.sph_sparta_atmospheric_parameters()
        if config['misc_query_databases']:
            self.sph_query_databases(timeout=config['misc_query_timeout'])
        self.sph_sparta_plot()


    def clean(self):
        '''
        Clean the reduction directory
        '''

        self._logger.info('====> Clean-up <====')

    
    def full_reduction(self):
        '''
        Performs a full reduction of a SPARTA data set
        '''
        
        self._logger.info('====> Full reduction <====')

        self.init_reduction()
        self.process_science()
        self.clean()
        
    ##################################################
    # SPHERE/SPARTA methods
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
        path = self.path

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
        file = open(Path(sphere.__file__).parent / 'instruments' / 'keywords_sparta.dat', 'r')
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

        # artificially add arm keyword
        files_info.insert(files_info.columns.get_loc('DPR TECH')+1, 'SEQ ARM', 'SPARTA')
            
        # drop files that are not handled, based on DPR keywords
        self._logger.debug('> drop unsupported file types')
        files_info.dropna(subset=['DPR TYPE'], inplace=True)
        files_info = files_info[(files_info['DPR TYPE'] == 'OBJECT,AO') & (files_info['OBS PROG ID'] != 'Maintenance')]

        # processed column
        files_info.insert(len(files_info.columns), 'PROCESSED', False)

        # convert times
        self._logger.debug('> convert times')
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=False)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=False)

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


    def sph_sparta_dtts(self, plot=True):
        '''
        Process SPARTA files for DTTS images

        Parameters
        ----------
        plot : bool
            Display and save diagnostic plot for quality check. Default is True
        '''
        
        self._logger.info('Process DTTS images')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_dtts', 
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info = self.files_info

        # build indices
        files = []
        img   = []
        for file, finfo in files_info.iterrows():
            self._logger.debug(f' * {file}')
            hdu = fits.open(f'{path.raw}/{file}.fits')
            
            ext  = hdu['IRPixelAvgFrame']
            NDIT = ext.header['NAXIS2']
            
            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))

            hdu.close()

        # create new dataframe
        self._logger.debug('> create data frame')
        dtts_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        dtts_info = dtts_info.align(files_info, level=0)[1]

        # extract data cube
        dtts_cube = np.zeros((len(dtts_info), 32, 32))
        nimg = 0
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')

            ext    = hdu['IRPixelAvgFrame']
            NDIT   = ext.header['NAXIS2']
            pixels = ext.data['Pixels'].reshape((-1, 32, 32))

            if NDIT:
                # timestamps
                time = Time(ext.data['Sec'] + ext.data['USec']*1e-6, format='unix')
                time.format = 'isot'
                dtts_info.loc[file, 'TIME']        = [str(t) for t in time]

                # DTTS images
                dtts_cube[nimg:nimg+NDIT] = pixels

                nimg += NDIT

            hdu.close()
            
        # updates times and compute timestamps
        toolbox.compute_times(dtts_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(dtts_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sph_sparta_dtts', sphere.ERROR)
            self._status = sphere.FATAL
            return
        
        # save
        self._logger.debug('> save dtts_frames.csv')
        dtts_info.to_csv(path.products / 'dtts_frames.csv')
        fits.writeto(path.products / 'dtts_cube.fits', dtts_cube, overwrite=True)
        self._dtts_info = dtts_info
        
        # plot
        if plot:
            self._logger.debug('> plot DTTS images')
            
            ncol  = 10
            nrow  = 10
            npage = int(np.ceil(nimg / (ncol*nrow)))+1
            vmax  = dtts_cube.max(axis=(1, 2)).mean()

            with PdfPages(path.products / 'dtts_images.pdf') as pdf:
                for page in range(npage):
                    self._logger.debug(f'  * page {page+1}/{npage}')

                    plt.figure(figsize=(3*ncol, 3*nrow))
                    plt.subplot(111)
                    
                    # master image
                    dtts_master = np.full((nrow*32, ncol*32), np.nan)
                    for row in range(nrow):
                        for col in range(ncol):
                            idx = page*nrow*ncol + row*ncol + col
                            
                            if idx < nimg:
                                xmin = col*32
                                xmax = (col+1)*32
                                ymin = (nrow-row-1)*32
                                ymax = (nrow-row)*32
                                
                                dtts_master[ymin:ymax, xmin:xmax] = dtts_cube[idx]

                                ts  = dtts_info['TIME'].values[idx]
                                date = ts[:10]
                                time = ts[11:]
                                plt.text(xmin+1, ymax-2, f'Date: {date}', size=14, weight='bold', color='w', ha='left', va='top', zorder=100)
                                plt.text(xmin+1, ymax-5, f'Time: {time}', size=14, weight='bold', color='w', ha='left', va='top', zorder=100)

                    plt.imshow(dtts_master, interpolation='nearest', vmin=0, vmax=vmax, cmap='inferno', zorder=0)

                    plt.xticks([])
                    plt.yticks([])

                    plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
                    
                    pdf.savefig()
                    plt.close()

        # update recipe execution
        self._update_recipe_status('sph_sparta_dtts', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE

        
    def sph_sparta_wfs_parameters(self):
        '''
        Process SPARTA files for Vis and IR WFS fluxes
        '''

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_wfs_parameters', 
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info = self.files_info

        #
        # VisLoop
        #
        
        self._logger.info('Process visible loop parameters')
        
        # build indices
        files = []
        img   = []
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')
            
            data = hdu['VisLoopParams']
            NDIT = data.header['NAXIS2']

            self._logger.debug(f' * {file} ==> {NDIT} records')

            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))

            hdu.close()

        # create new dataframe
        self._logger.debug('> create data frame')
        visloop_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        visloop_info = visloop_info.align(files_info, level=0)[1]

        # extract data
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')

            ext  = hdu['VisLoopParams']
            NDIT = ext.header['NAXIS2']
            
            if NDIT:
                # timestamps
                time = Time(ext.data['Sec'] + ext.data['USec']*1e-6, format='unix')
                time.format = 'isot'
                visloop_info.loc[file, 'TIME'] = [str(t) for t in time]

                # VisLoop parameters
                visloop_info.loc[file, 'focus_avg']      = ext.data['Focus_avg']
                visloop_info.loc[file, 'TTx_avg']        = ext.data['TTx_avg']
                visloop_info.loc[file, 'TTy_avg']        = ext.data['TTy_avg']
                visloop_info.loc[file, 'DMPos_avg']      = ext.data['DMPos_avg']
                visloop_info.loc[file, 'ITTMPos_avg']    = ext.data['ITTMPos_avg']
                visloop_info.loc[file, 'DMSatur_avg']    = ext.data['DMSatur_avg']
                visloop_info.loc[file, 'DMAberr_avg']    = ext.data['DMAberr_avg']
                visloop_info.loc[file, 'flux_total_avg'] = ext.data['Flux_avg']
                
            hdu.close()

        # convert VisWFS flux in photons per subaperture. Flux_avg is the flux on the whole pupil made of 1240 subapertures
        photon_to_ADU = 17     # from Jean-François Sauvage
        gain = visloop_info['AOS VISWFS MODE'].str.split('_', expand=True)[1].astype(float)
        visloop_info['flux_subap_avg'] = visloop_info['flux_total_avg'] / gain / photon_to_ADU  # in photons per subaperture

        # updates times and compute timestamps
        toolbox.compute_times(visloop_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(visloop_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sph_sparta_wfs_parameters', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # save
        self._logger.debug('> save visloop_info.csv')
        visloop_info.to_csv(path.products / 'visloop_info.csv')
        self._visloop_info = visloop_info
    
        #
        # IRLoop
        #

        self._logger.info('Process IR loop parameters')
        
        # build indices
        files = []
        img   = []
        for file, finfo in files_info.iterrows():
            self._logger.debug(f' * {file}')
            hdu = fits.open(f'{path.raw}/{file}.fits')
            
            data = hdu['IRLoopParams']
            NDIT = data.header['NAXIS2']
            
            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))

            hdu.close()

        # create new dataframe
        self._logger.debug('> create data frame')
        irloop_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        irloop_info = irloop_info.align(files_info, level=0)[1]

        # extract data
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')

            ext  = hdu['IRLoopParams']
            NDIT = ext.header['NAXIS2']

            if NDIT:
                # timestamps
                time = Time(ext.data['Sec'] + ext.data['USec']*1e-6, format='unix')
                time.format = 'isot'
                irloop_info.loc[file, 'TIME'] = [str(t) for t in time]
                
                # VisLoop parameters
                irloop_info.loc[file, 'DTTPPos_avg'] = ext.data['DTTPPos_avg']
                irloop_info.loc[file, 'DTTPRes_avg'] = ext.data['DTTPRes_avg']
                irloop_info.loc[file, 'flux_avg']    = ext.data['Flux_avg']
                
            hdu.close()
    
        # updates times and compute timestamps
        toolbox.compute_times(irloop_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(irloop_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sph_sparta_wfs_parameters', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # save
        self._logger.debug('> save irloop_info.csv')
        irloop_info.to_csv(path.products / 'irloop_info.csv')
        self._irloop_info = irloop_info

        # update recipe execution
        self._update_recipe_status('sph_sparta_wfs_parameters', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
        

    def sph_sparta_atmospheric_parameters(self):
        '''
        Process SPARTA files for atmospheric parameters
        '''

        self._logger.info('Process atmospheric parameters')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_atmospheric_parameters', 
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info = self.files_info
    
        #
        # Atmospheric parameters
        #
        
        # build indices
        files = []
        img   = []
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')
            
            data = hdu['AtmPerfParams']
            NDIT = data.header['NAXIS2']

            self._logger.debug(f' * {file} ==> {NDIT} records')

            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))

            hdu.close()

        # create new dataframe
        self._logger.debug('> create data frame')
        atmos_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        atmos_info = atmos_info.align(files_info, level=0)[1]

        # extract data
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')

            ext  = hdu['AtmPerfParams']
            NDIT = ext.header['NAXIS2']
            
            if NDIT:
                # timestamps
                time = Time(ext.data['Sec'] + ext.data['USec']*1e-6, format='unix')
                time.format = 'isot'
                atmos_info.loc[file, 'TIME'] = [str(t) for t in time]

                # Atm parameters
                atmos_info.loc[file, 'r0']        = ext.data['R0']
                atmos_info.loc[file, 'windspeed'] = ext.data['WindSpeed']
                atmos_info.loc[file, 'strehl']    = ext.data['StrehlRatio']
                
            hdu.close()

        # updates times and compute timestamps
        toolbox.compute_times(atmos_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(atmos_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sph_sparta_atmospheric_parameters', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # remove bad values
        atmos_info.loc[np.logical_or(atmos_info['strehl'] <= 0.05, atmos_info['strehl'] > 0.98), 'strehl'] = np.nan
        atmos_info.loc[np.logical_or(atmos_info['r0'] <= 0, atmos_info['r0'] > 0.9), 'r0'] = np.nan
        atmos_info.loc[np.logical_or(atmos_info['windspeed'] <= 0, atmos_info['windspeed'] > 50), 'windspeed'] = np.nan

        # tau0 and the seeing from r0
        atmos_info['tau0']   = 0.314*atmos_info['r0'] / atmos_info['windspeed']
        atmos_info['seeing'] = np.rad2deg(wave_wfs / atmos_info['r0']) * 3600

        # IMPLEMENT:
        # we compute the zenith seeing: seeing(zenith) = seeing(AM)/AM^(3/5)
        atmos_info['seeing_zenith'] = atmos_info['seeing'] / np.power(atmos_info['AIRMASS'], 3/5)
        atmos_info['r0_zenith']     = atmos_info['r0'] * np.power(atmos_info['AIRMASS'], 3/5)
        atmos_info['tau0_zenith']   = atmos_info['tau0'] * np.power(atmos_info['AIRMASS'], 3/5)
        
        # save
        self._logger.debug('> save atmospheric_info.csv')
        atmos_info.to_csv(path.products / 'atmospheric_info.csv')
        self._atmos_info = atmos_info
        
        # update recipe execution
        self._update_recipe_status('sph_sparta_atmospheric_parameters', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
        

    def sph_query_databases(self, timeout=5):
        '''
        Query ESO databases for additional atmospheric information

        See details on the ESO webpage:
            http://archive.eso.org/cms/eso-data/ambient-conditions/paranal-ambient-query-forms.html

        The following instruments are queried:
            - DIMM: Differential Image Moption Monitor
            - MASS: Multi-Aperture Scintillation Sensor
            - SLODAR: SLOpe Detection And Ranging
            - LHATPRO: Low Humidity And Temperature PROfiling microwave radiometer
            - METEO

        Parameters
        ----------
        timeout : float
            Network request timeout, in seconds. Default is 5
        '''
        
        self._logger.info('Query ESO databases')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_query_databases', 
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info = self.files_info

        # times
        time_start = Time(files_info['DATE'].min()).isot
        time_end   = Time(files_info['DATE-OBS'].max()).isot

        #
        # MASS-DIMM
        #
        self._logger.info('Query MASS-DIMM')

        url = f'http://archive.eso.org/wdb/wdb/asm/mass_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_fwhm=1&tab_fwhmerr=0&tab_tau=1&tab_tauerr=0&tab_tet=0&tab_teterr=0&tab_alt=0&tab_alterr=0&tab_fracgl=1&tab_turbfwhm=1&tab_tau0=1&tab_tet0=0&tab_turb_alt=0&tab_turb_speed=1'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                mass_dimm_info = pd.read_csv(data, index_col=0, comment='#')
                mass_dimm_info.rename(columns={'Date time': 'date',
                                               'MASS Tau0 [s]': 'mass_tau0',
                                               'MASS-DIMM Cn2 fraction at ground': 'mass-dimm_Cn2_frac_gl',
                                               'MASS-DIMM Tau0 [s]': 'mass-dimm_tau0',
                                               'MASS-DIMM Turb Velocity [m/s]': 'mass-dimm_turb_speed',
                                               'MASS-DIMM Seeing ["]': 'mass-dimm_seeing',
                                               'Free Atmosphere Seeing ["]': 'mass_freeatmos_seeing'}, inplace=True)
                mass_dimm_info.to_csv(path.products / 'mass-dimm_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to MASS-DIMM timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty MASS-DIMM data')

        #
        # DIMM
        #
        self._logger.info('Query DIMM')

        if Time(time_end) <= Time('2016-04-05T00:00:00.000'):
            url = f'http://archive.eso.org/wdb/wdb/asm/historical_ambient_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_fwhm=1&tab_airmass=0&tab_rfl=0&tab_tau=1&tab_tet=0&top=1000000'
        else:
            url = f'http://archive.eso.org/wdb/wdb/asm/dimm_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_fwhm=1&tab_rfl=0&tab_rfl_time=0'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                dimm_info = pd.read_csv(data, index_col=0, comment='#')
                dimm_info.rename(columns={'Date time': 'date',
                                          'DIMM Seeing ["]': 'dimm_seeing',
                                          'Tau0 [s]': 'old_dimm_tau0'}, inplace=True)
                dimm_info.to_csv(path.products / 'dimm_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to DIMM timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty DIMM data')

        #
        # SLODAR
        #

        self._logger.info('Query SLODAR')

        url = f'http://archive.eso.org/wdb/wdb/asm/slodar_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_cnsqs_uts=1&tab_fracgl300=1&tab_fracgl500=1&tab_hrsfit=1&tab_fwhm=1'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                slodar_info = pd.read_csv(data, index_col=0, comment='#')
                slodar_info.rename(columns={'Date time': 'date', 'Cn2 above UTs [10**(-15)m**(1/3)]': 'Cn2_above_UT',
                                            'Cn2 fraction below 300m': 'slodar_frac_below_300m',
                                            'Cn2 fraction below 500m': 'slodar_frac_below_500m',
                                            'Surface layer profile [10**(-15)m**(1/3)]': 'slodar_surface_layer',
                                            'Seeing ["]': 'slodar_seeing'}, inplace=True)

                wave_num = 2*np.pi/wave_wfs
                slodar_info['slodar_r0_above_UT'] = np.power(0.423*(wave_num**2)*slodar_info['Cn2_above_UT']*1.e-15, -3/5)
                slodar_info['slodar_seeing_above_UT'] = np.rad2deg(wave_wfs/slodar_info['slodar_r0_above_UT'])*3600
                slodar_info['slodar_Cn2_total'] = np.power(slodar_info['slodar_seeing']/2e7, 1/0.6)  # in m^1/3
                slodar_info['slodar_surface_layer_frac'] = slodar_info['slodar_surface_layer']*1e-15 / slodar_info['slodar_Cn2_total']
                slodar_info.to_csv(path.products / 'slodar_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to SLODAR timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty SLODAR data')

        #
        # Telescope seeing
        #

        self._logger.info('Query telescope seeing from SPHERE data')

        url = f'http://archive.eso.org/wdb/wdb/eso/sphere/query?wdbo=csv&night={time_start:s}..{time_end:s}&tab_prog_id=0&tab_dp_id=0&tab_ob_id=0&tab_exptime=0&tab_dp_cat=0&tab_tpl_start=0&tab_dp_type=0&tab_dp_tech=0&tab_seq_arm=0&tab_ins3_opti5_name=0&tab_ins3_opti6_name=0&tab_ins_comb_vcor=0&tab_ins_comb_iflt=0&tab_ins_comb_pola=0&tab_ins_comb_icor=0&tab_det_dit1=0&tab_det_seq1_dit=0&tab_det_ndit=0&tab_det_read_curname=0&tab_ins2_opti2_name=0&tab_det_chip_index=0&tab_ins4_comb_rot=0&tab_ins1_filt_name=0&tab_ins1_opti1_name=0&tab_ins1_opti2_name=0&tab_ins4_opti11_name=0&tab_tel_ia_fwhm=1&tab_tel_ia_fwhmlin=1&tab_tel_ia_fwhmlinobs=1&tab_tel_ambi_windsp=0&tab_night=1&tab_fwhm_avg=0&top=1000'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                sphere_info = pd.read_csv(data, index_col=6, comment='#')
                
                keys_to_drop = ['Release Date', 'Object', 'RA', 'DEC', 'Target Ra Dec', 'Target l b', 'OBS Target Name']
                for key in keys_to_drop:
                    sphere_info.drop(key, axis=1, inplace=True)
                    
                sphere_info.to_csv(path.products / 'ambi_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to SPHERE archive timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty SPHERE archive data')

        #
        # Meteo tower for wind speed/direction and temperature
        #

        self._logger.info('Query meteo tower')

        url = f'http://archive.eso.org/wdb/wdb/asm/meteo_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_press=0&tab_presqnh=0&tab_temp1=1&tab_temp2=0&tab_temp3=0&tab_temp4=0&tab_tempdew1=0&tab_tempdew2=0&tab_tempdew4=0&tab_dustl1=0&tab_dustl2=0&tab_dusts1=0&tab_dusts2=0&tab_rain=0&tab_rhum1=0&tab_rhum2=0&tab_rhum4=0&tab_wind_dir1=1&tab_wind_dir1_180=0&tab_wind_dir2=0&tab_wind_dir2_180=0&tab_wind_speed1=1&tab_wind_speed2=0&tab_wind_speedu=0&tab_wind_speedv=0&tab_wind_speedw=0'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                meteo_info = pd.read_csv(data, index_col=0, comment='#')
                meteo_info.rename(columns={'Date time': 'date',
                                           'Air Temperature at 30m [C]': 'air_temperature_30m',
                                           'Wind Direction at 30m (0/360) [deg]': 'winddir_30m',
                                           'Wind Speed at 30m [m/s]': 'windspeed_30m'}, inplace=True)
                meteo_info.to_csv(path.products / 'meteo_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to meteo tower timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty meteo tower data')
        
        #
        # LATHPRO
        #

        self._logger.info('Query LATHPRO')

        url = f'http://archive.eso.org/wdb/wdb/asm/lhatpro_paranal/query?wdbo=csv&start_date={time_start:s}..{time_end:s}&tab_integration=0&tab_lwp0=0'

        try:
            req = requests.get(url, timeout=timeout)
            if req.status_code == requests.codes.ok:
                self._logger.debug('  ==> Valid response received')

                data = io.StringIO(req.text)
                lathpro_info = pd.read_csv(data, index_col=0, comment='#')
                lathpro_info.rename(columns={'Date time': 'date',
                                             'IR temperature [Celsius]': 'lathpro_IR_temperature',
                                             'Precipitable Water Vapour [mm]': 'lathpro_pwv'}, inplace=True)
                lathpro_info.to_csv(path.products / 'lathpro_info.csv')
            else:
                self._logger.debug('  ==> Query error')
        except requests.ReadTimeout:
            self._logger.error('  ==> Request to meteo tower timed out')
        except pd.errors.EmptyDataError:
            self._logger.error('  ==> Empty meteo tower data')

        # update recipe execution
        self._update_recipe_status('sph_query_databases', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
    

    def sph_sparta_plot(self):
        '''
        Plot results of the SPARTA file analysis
        '''

        self._logger.info('Plot SPARTA and atmospheric data')
        
        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_plot',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path
        files_info   = self.files_info
        visloop_info = self.visloop_info
        irloop_info  = self.irloop_info
        atmos_info   = self.atmospheric_info

        # make sure we have times where needed
        self._logger.debug('> convert times')
        atmos_info['TIME']   = pd.to_datetime(atmos_info['TIME'], utc=False)
        visloop_info['TIME'] = pd.to_datetime(visloop_info['TIME'], utc=False)
        irloop_info['TIME']  = pd.to_datetime(irloop_info['TIME'], utc=False)
        
        # additional databases
        file = path.products / 'mass-dimm_info.csv'
        if file.exists():
            mass_dimm_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            mass_dimm_info = None

        file = path.products / 'dimm_info.csv'
        if file.exists():
            dimm_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            dimm_info = None

        file = path.products / 'slodar_info.csv'
        if file.exists():
            slodar_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            slodar_info = None

        file = path.products / 'ambi_info.csv'
        if file.exists():
            sphere_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            sphere_info = None

        file = path.products / 'meteo_info.csv'
        if file.exists():
            meteo_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            meteo_info = None

        file = path.products / 'lathpro_info.csv'
        if file.exists():
            lathpro_info = pd.read_csv(file, index_col=0, parse_dates=True)
        else:
            lathpro_info = None

        # times
        time_start = Time(files_info['DATE'].min())
        time_end   = Time(files_info['DATE-OBS'].max())

        #
        # plot
        #
        dateFormatter = mdates.DateFormatter('%H:%M')
        
        fig = plt.figure(1, figsize=(12, 10))
        plt.rcParams.update({'font.size': 14})

        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1], width_ratios=[1, 1])

        # seeing
        ax = plt.subplot(gs[0, 0])

        ax.plot_date(Time(atmos_info['TIME']).plot_date, atmos_info['seeing_zenith'], '.', color='darkorange', markeredgecolor='none', label='SPARTA')
        if mass_dimm_info is not None:
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass-dimm_seeing'], '.', color='palevioletred', markeredgecolor='none', label='MASS-DIMM')
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass_freeatmos_seeing'], '.', color='lime', markeredgecolor='none', label='MASS')
        if dimm_info is not None:
            ax.plot_date(Time(dimm_info.index).plot_date, dimm_info['dimm_seeing'], '.', color='dimgrey', markeredgecolor='none', label='DIMM')
        if slodar_info is not None:
            ax.plot_date(Time(slodar_info.index).plot_date, slodar_info['slodar_seeing_above_UT'], '.', color='magenta', markeredgecolor='none', label='SLODAR above UT')
        if sphere_info is not None:
            ax.plot_date(Time(sphere_info.index).plot_date, sphere_info['TEL IA FWHMLIN'], '.', color='rosybrown', markeredgecolor='none', label='TEL.IA.FWHMLIN')

        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylabel('Seeing ["]')
        ax.xaxis.set_ticklabels([])
        # ax.xaxis.set_major_formatter(dateFormatter)
        
        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')
        
        # r0
        ax = plt.subplot(gs[0, 1])

        ax.plot_date(Time(atmos_info['TIME']).plot_date, atmos_info['tau0_zenith']*1000, '.', color='darkgreen', markeredgecolor='none', label='SPARTA')
        if mass_dimm_info is not None:
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass-dimm_tau0']*1000, '.', color='palevioletred', label='MASS-DIMM', markeredgecolor='none')
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass_tau0']*1000., '.', color='dimgrey', label='MASS', markeredgecolor='none')

        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylabel('$\\tau_0$ [ms]')
        ax.xaxis.set_ticklabels([])
        # ax.xaxis.set_major_formatter(dateFormatter)

        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')
            
        # Vis and NIR WFS fluxes
        ax = plt.subplot(gs[1, 0])

        flux_visloop = visloop_info['flux_subap_avg']
        flux_irloop  = irloop_info['flux_avg']
        
        ax.plot_date(Time(visloop_info['TIME']).plot_date, flux_visloop, '.', color='blue', markeredgecolor='none', label='Vis WFS')
        ax.plot_date(Time(irloop_info['TIME']).plot_date, flux_irloop, '.', color='red', markeredgecolor='none', label='IR DTTS')

        ymin = np.percentile(np.append(flux_irloop, flux_visloop), 10)/10
        ymax = np.percentile(np.append(flux_irloop, flux_visloop), 90)*10
        
        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylim(ymin, ymax)
        ax.set_yscale('log', nonpositive='clip')
        ax.set_ylabel('Flux in photons/aperture')
        ax.xaxis.set_ticklabels([])
        # ax.xaxis.set_major_formatter(dateFormatter)

        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')
        
        # Strehl
        ax = plt.subplot(gs[1, 1])

        ax.plot_date(Time(atmos_info['TIME']).plot_date, atmos_info['strehl']*100, '.', color='darkorchid', markeredgecolor='none', label='SPARTA')

        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylim(0, 100)
        ax.set_ylabel('Strehl ratio [%]')
        ax.xaxis.set_ticklabels([])
        # ax.xaxis.set_major_formatter(dateFormatter)

        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')
        
        # ground-layer fraction
        ax = plt.subplot(gs[2, 0])

        if mass_dimm_info is not None:
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass-dimm_Cn2_frac_gl']*100, '.', color='palevioletred', label='MASS-DIMM', markeredgecolor='none')
        if slodar_info is not None:
            ax.plot_date(Time(slodar_info.index).plot_date, slodar_info['slodar_surface_layer_frac']*100, '.', color='black', markeredgecolor='none', label='SLODAR surface layer')
            ax.plot_date(Time(slodar_info.index).plot_date, slodar_info['slodar_frac_below_500m']*100, '.', color='red', markeredgecolor='none', label='SLODAR <500 m')
            ax.plot_date(Time(slodar_info.index).plot_date, slodar_info['slodar_frac_below_300m']*100, '.', color='blue', markeredgecolor='none', label='SLODAR <300 m')

        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylim(0, 100)
        ax.set_ylabel('Ground layer fraction [%]')
        ax.xaxis.set_major_formatter(dateFormatter)

        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')

        # wind speed
        ax = plt.subplot(gs[2, 1])

        ax.plot_date(Time(atmos_info['TIME']).plot_date, atmos_info['windspeed'], '.', color='darkgreen', markeredgecolor='none', label='SPARTA')
        if mass_dimm_info is not None:
            ax.plot_date(Time(mass_dimm_info.index).plot_date, mass_dimm_info['mass-dimm_turb_speed'], '.', color='palevioletred', label='MASS', markeredgecolor='none')
        if meteo_info is not None:
            ax.plot_date(Time(meteo_info.index).plot_date, meteo_info['windspeed_30m'], '.', color='rosybrown', label='ASM 30 m', markeredgecolor='none')

        ax.tick_params(axis='x', labelrotation=45, labelsize='small')
        ax.set_xlim(time_start.plot_date, time_end.plot_date)
        ax.set_ylabel('Wind speed [m/s]')
        ax.xaxis.set_major_formatter(dateFormatter)
        
        ax.grid(True)
        ax.legend(frameon=True, loc='best', fontsize='x-small')

        # tweak
        plt.subplots_adjust(left=0.07, right=0.98, bottom=0.07, top=0.98, wspace=0.25, hspace=0.12)

        plt.savefig(path.products / 'sparta_plot.pdf')
        
        # update recipe execution
        self._update_recipe_status('sph_query_databases', sphere.SUCCESS)

        # reduction status
        self._status = sphere.COMPLETE
    

    def sph_sparta_clean(self, delete_raw=False, delete_products=False):
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
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_clean',
                                         self.recipe_requirements, logger=self._logger):
            return

        # parameters
        path = self.path

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
        self._update_recipe_status('sph_sparta_clean', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
