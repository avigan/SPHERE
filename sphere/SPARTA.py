import pandas as pd
import logging
import numpy as np
import collections
import configparser
import shutil
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.time import Time
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages

import sphere
import sphere.utils as utils
import sphere.toolbox as toolbox

_log = logging.getLogger(__name__)


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
    recipe_requirements = {
        'sort_files': [],
        'sph_sparta_dtts': ['sort_files'],
        'sph_sparta_wfs_parameters': ['sort_files'],
        'sph_sparta_atmospheric_parameters': ['sort_files'],
        'sph_sparta_query_databases': ['sort_file', 'sph_sparta_dtts'],
        'sph_ifs_clean': []
    }

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
        # reduction status
        #
        reduction._status = sphere.INIT
        reduction._recipes_status = collections.OrderedDict()
        
        # reload any existing data frames
        # FIXME: to be implemented
        # reduction._read_info()

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
    def recipe_status(self):
        return self._recipes_status

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

            # update instrument mode
            self._mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 MODE'][0]
        else:
            files_info = None

        # save data frames in instance variables
        self._files_info = files_info

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
        self.sph_spart_wfs_flux()
        self.sph_sparta_atmospheric_parameters()

        if config['misc_query_database']:
            self.sph_sparta_query_databases(timeout=config['misc_query_timeout'])

        
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

        # artificially add arm and telescope keywords
        files_info.insert(files_info.columns.get_loc('DPR TECH')+1, 'SEQ ARM', 'SPARTA')
        files_info.insert(files_info.columns.get_loc('INS4 DROT2 RA'), 'TEL GEOLON', -26.6268)
        files_info.insert(files_info.columns.get_loc('INS4 DROT2 RA'), 'TEL GEOLAT', -70.4045)
        files_info.insert(files_info.columns.get_loc('INS4 DROT2 RA'), 'TEL GEOELEV', 2648)
            
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
        self._logger.debug('> create frames_info data frame')
        dtts_frames_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        dtts_frames_info = dtts_frames_info.align(files_info, level=0)[1]

        # extract data cube
        dtts_cube = np.zeros((len(dtts_frames_info), 32, 32))
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
                dtts_frames_info.loc[file, 'TIME']        = [str(t) for t in time]

                # DTTS images
                dtts_cube[nimg:nimg+NDIT] = pixels

                nimg += NDIT

            hdu.close()
            
        # updates times and compute timestamps
        toolbox.compute_times(dtts_frames_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(dtts_frames_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sort_frames', sphere.ERROR)
            self._status = sphere.FATAL
            return
        
        # save
        self._logger.debug('> save dtts_frames.csv')
        dtts_frames_info.to_csv(path.products / 'dtts_frames.csv')
        fits.writeto(path.products / 'dtts_cube.fits', dtts_cube, overwrite=True)

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

                                ts  = dtts_frames_info['TIME'].values[idx]
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
        self._update_recipe_status('sph_sparta_process', sphere.SUCCESS)

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
        self._logger.debug('> create frames_info data frame')
        visloop_frames_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        visloop_frames_info = visloop_frames_info.align(files_info, level=0)[1]

        # extract data
        for file, finfo in files_info.iterrows():
            hdu = fits.open(f'{path.raw}/{file}.fits')

            ext  = hdu['VisLoopParams']
            NDIT = ext.header['NAXIS2']
            
            if NDIT:
                # timestamps
                time = Time(ext.data['Sec'] + ext.data['USec']*1e-6, format='unix')
                time.format = 'isot'
                visloop_frames_info.loc[file, 'TIME']        = [str(t) for t in time]

                # VisLoop parameters
                visloop_frames_info.loc[file, 'Focus_avg']   = ext.data['Focus_avg']
                visloop_frames_info.loc[file, 'TTx_avg']     = ext.data['TTx_avg']
                visloop_frames_info.loc[file, 'TTy_avg']     = ext.data['TTy_avg']
                visloop_frames_info.loc[file, 'DMPos_avg']   = ext.data['DMPos_avg']
                visloop_frames_info.loc[file, 'ITTMPos_avg'] = ext.data['ITTMPos_avg']
                visloop_frames_info.loc[file, 'DMSatur_avg'] = ext.data['DMSatur_avg']
                visloop_frames_info.loc[file, 'DMAberr_avg'] = ext.data['DMAberr_avg']
                
            hdu.close()
        
        # updates times and compute timestamps
        toolbox.compute_times(visloop_frames_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(visloop_frames_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sort_frames', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # save
        self._logger.debug('> save visloop_frames.csv')
        visloop_frames_info.to_csv(path.products / 'visloop_frames.csv')

        stop
        
        #
        # IRLoop
        #

        self._logger.info('Process IR loop parameters')
        
        # build indices
        files = []
        img   = []
        times = []
        nimg  = 0
        for file, finfo in files_info.iterrows():
            self._logger.debug(f' * {file}')
            hdu = fits.open(f'{path.raw}/{file}.fits')
            
            data = hdu['IRLoopParams']
            NDIT = data.header['NAXIS2']
            time = Time(data.data['Sec'] + data.data['USec']*1e-6, format='unix')
            time.format = 'isot'
            nimg += NDIT
            
            files.extend(np.repeat(file, NDIT))
            img.extend(list(np.arange(NDIT)))
            times.extend([str(t) for t in time])

            hdu.close()

        # create new dataframe
        self._logger.debug('> create frames_info data frame')
        irloop_frames_info = pd.DataFrame(columns=files_info.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))

        # expand files_info into frames_info
        irloop_frames_info = irloop_frames_info.align(files_info, level=0)[1]

        # updates times and compute timestamps
        irloop_frames_info['TIME'] = times
        toolbox.compute_times(irloop_frames_info, logger=self._logger)

        # compute angles (ra, dec, parang)
        ret = toolbox.compute_angles(irloop_frames_info, logger=self._logger)
        if ret == sphere.ERROR:
            self._update_recipe_status('sort_frames', sphere.ERROR)
            self._status = sphere.FATAL
            return

        # save
        self._logger.debug('> save irloop_frames.csv')
        irloop_frames_info.to_csv(path.products / 'irloop_frames.csv')

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
        # IMPLEMENTATION
        #

        # update recipe execution
        self._update_recipe_status('sph_sparta_atmospheric_parameters', sphere.SUCCESS)

        # reduction status
        self._status = sphere.INCOMPLETE
        

    def sph_sparta_query_databases(self, timeout=5):
        '''
        Query ESO databases for additional atmospheric information

        Parameters
        ----------
        timeout : float
            Network request timeout, in seconds. Default is 5
        '''
        
        self._logger.info('Query ESO databases')

        # check if recipe can be executed
        if not toolbox.recipe_executable(self._recipes_status, self._status, 'sph_sparta_query_databases', 
                                         self.recipe_requirements, logger=self._logger):
            return

        #
        # TO BE IMPLEMENTED
        #

        # update recipe execution
        self._update_recipe_status('sph_sparta_query_databases', sphere.SUCCESS)

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
