import os
import glob
import shutil
import math
import logging
import numpy as np
import pandas as pd
import xml.etree.ElementTree as etree

import sphere.IRDIS as IRDIS
import sphere.IFS as IFS

from pathlib import Path
from astropy.io import fits
from astropy.time import Time

_log = logging.getLogger(__name__)


def process_mainFiles(mainFiles, files, logger=_log):
    '''
    Process top-level file association XML from the ESO archive

    Parameters
    ----------
    mainFiles : etree element
        Main file element

    files : list
        List where the files will be appended

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''
    # append file names to the list
    for file in mainFiles:
        fname = file.attrib['name']
        files.append(fname)

        logger.debug(' ==> {0}'.format(fname))
        

def process_association(tree, files, logger=_log):
    '''
    Process file association XML from the ESO archive

    Parameters
    ----------
    tree : etree element
        Associated file element

    files : list
        List where the files will be appended

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''
    catg = tree.attrib['category']

    logger.debug(catg)

    # skip unused calibrations
    if (catg == 'IFS_STD_ASTROM') or (catg == 'IFS_STD_PHOT') or \
       (catg == 'IFS_DIST') or (catg == 'IRD_CLI_PHOT') or \
       (catg == 'IRD_DIST'):
        logger.debug(' ==> skipping')
        return

    # process differently mainFiles from associatedFiles
    for elt in tree:
        if elt.tag == 'mainFiles':
            logger.debug('mainFiles')
            process_mainFiles(elt, files)
        elif elt.tag == 'associatedFiles':
            logger.debug('associatedFiles')
            for nelt in elt:
                process_association(nelt, files, logger=logger)


def sort_files_from_xml(path, logger=_log):
    '''Sort files downloaded from the ESO archive with associated raw
       calibrations

    When choosing this option from the archive, the data package
    includes xml files that provide the association between the
    science files and the calibration files. The method uses these
    xml files to copy the science and associated raw files into
    dedicated directories.

    For a given executed OB, the method will extract the target
    name from the OBS.NAME keyword and the OB ID from the OBS.ID
    keyword. It will create a directory {OBS.NAME}_id={OBS.ID} and
    copy the FITS files. Note that any space in the name will be
    replaced by an underscore to avoid subsequent issues with
    esorex in the data reduction.

    The files are *copied* because science data acquired on the
    same night with the same setup can share some calibrations. It
    is therefore safer to copy the needed fiels instead of moving
    them. At the end, all the original files are more to the
    all_files/ subdirectory.

    Parameters
    ----------
    path : str
        Path where to look for XML files

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''

    xml_files = list(path.glob('*.xml'))

    logger.info('Sort data based on XML files (ESO automated calibration selection)')
    logger.info(' * {0} XML files\n'.format(len(xml_files)))
    
    # sort files
    for file in xml_files:
        tree = etree.parse(file)
        root = tree.getroot()            

        logger.info(' * {}'.format(file.name))
        
        # process only IFS and IRDIS science data
        catg = root.attrib['category']
        if (catg.find('ACQUISITION') != -1):
            continue

        # get target name from first mainFile element
        scifiles = root.find('mainFiles')
        filename = scifiles[0].attrib['name']
        
        # Mac OS X replaces : by _ in file names...
        if not (path / '{}.fits'.format(filename)).exists():
            filename = filename.replace(':', '_')
        
        if not (path / '{}.fits'.format(filename)).exists():
            logger.info('   ==> file {} does not exist. Skipping'.format(filename))
            continue

        fpath = path / '{}.fits'.format(filename)
        hdr = fits.getheader(fpath)
        
        # target and arm
        target = hdr['HIERARCH ESO OBS NAME']
        obs_id = hdr['HIERARCH ESO OBS ID']
        mjd    = Time(math.floor(float(hdr['MJD-OBS']+0.5))-0.5, format='mjd')
        night  = mjd.isot[:10]
        if catg == 'SCIENCE_OBJECT_AO':
            instrument = 'SPARTA'
        else:
            try:
                arm = hdr['HIERARCH ESO SEQ ARM']
            except KeyError:
                logger.error('No \'HIERARCH ESO SEQ ARM\' keyword in {}'.format(fpath))
                continue
                
            if arm == 'IRDIS':
                instrument = 'IRDIS'
            elif arm == 'IFS':
                instrument = 'IFS'                
            else:
                logger.error('Unknown arm {0}'.format(arm))
                continue

        # get files
        files = []
        process_association(root, files, logger=logger)

        # target path
        directory = '{0}_id={1}'.format(target, obs_id)
        directory = '_'.join(directory.split())
        target_path = path / directory / night / instrument / 'raw'
        target_path.mkdir(parents=True, exist_ok=True)

        # copy files
        for filename in files:
            fpath = path / '{}.fits'.format(filename)
            
            # Mac OS X replaces : by _ in file names...
            if not fpath.exists():
                filename = filename.replace(':', '_')
                fpath  = path / '{}.fits'.format(filename)
            
            # check if file actually exists
            if not fpath.exists():
                logger.info(' ==> file {} does not exist. Skipping.'.format(fpath))
                continue
            
            # copy if needed
            nfpath = target_path / '{}.fits'.format(filename)
            if not nfpath.exists():
                shutil.copy(fpath, nfpath)

        # print status
        logger.debug('{0} - id={1}'.format(target, obs_id))
        logger.debug(' ==> found {0} files'.format(len(files)))
        logger.debug(' ==> copied to {0}'.format(target_path))

    # move all files
    path_new = path / 'all_files'
    path_new.mkdir(parents=True, exist_ok=True)

    files = []
    files.extend(list(path.glob('*.fits')))
    files.extend(list(path.glob('*.xml')))
    files.extend(list(path.glob('*.txt')))

    if len(files) != 0:
        for file in files:
            file.rename(path_new / file.name)


def sort_files_from_fits(path, logger=_log):
    '''Sort FITS files based only based on their headers

    Contrary to sort_files_from_xml(), this method is dumb in the
    sense that it does not try to keep any file associated. It
    just sorts the FITS files as a function of the
    sub-system. Still very convenient when manually downloading
    data and calibrations.

    And contrary to sort_files_from_xml(), the files are directly
    moved to their sub-system directory. This makes this method
    much faster than sort_files_from_xml(). At the end, all
    unsorted files are moved to an unsorted_files/ subdirectory.

    Parameters
    ----------
    path : str
        Path where to look for FITS files

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    '''

    fits_files = list(path.glob('*.fits'))

    logger.info('Sort data based on FITS files')
    logger.info(' * {0} FITS files\n'.format(len(fits_files)))

    # sort files
    for file in fits_files:
        logger.info(' * {}'.format(file.name))
        
        # target and arm
        hdr = fits.getheader(file)

        try:
            target = hdr['HIERARCH ESO OBS NAME']
            obs_id = hdr['HIERARCH ESO OBS ID']
            dpr_type = hdr['HIERARCH ESO DPR TYPE']
        except KeyError:
            logger.error('Missing ESO HIERARCH keywords in {}'.format(file))
            continue

        if dpr_type == 'OBJECT,AO':
            instrument = 'SPARTA'
        else:
            try:
                arm = hdr['HIERARCH ESO SEQ ARM']
            except KeyError:
                logger.error('No \'HIERARCH ESO SEQ ARM\' keyword in {}'.format(file))
                continue
            
            if arm == 'IRDIS':
                instrument = 'IRDIS'
            elif arm == 'IFS':
                instrument = 'IFS'                
            else:
                logger.error('Unknown arm {0}'.format(arm))
                continue

        # target path
        target_path = path / instrument / 'raw'
        target_path.mkdir(parents=True, exist_ok=True)

        # move file
        file.rename(target_path / file.name)

        # print status
        logger.debug('{0} - id={1}'.format(target, obs_id))
        logger.debug(' ==> copied to {0}'.format(target_path))

    # move all files
    path_new = path / 'unsorted_files'
    path_new.mkdir(parents=True, exist_ok=True)

    files = []
    files.extend(list(path.glob('*.fits')))
    files.extend(list(path.glob('*.txt')))

    if len(files) != 0:
        for file in files:
            file.rename(path_new / file.name)

            
def classify_irdis_dataset(path, logger=_log):
    '''Classify an IRDIS dataset based on the science files

    Parameters
    ----------
    path : str
        Path to the directory containing the dataset

    logger : logHandler object
        Log handler for the reduction. Default is root logger

    Returns
    -------
    mode : str
        Generic string representing the name of the mode. None in case
        of failure.

    '''
    
    # zeroth-order reduction validation
    raw = path / 'raw'
    if not raw.exists():
        logger.error('No raw/ subdirectory. {0} is not a valid reduction path!'.format(path))
        return None

    # list all fits files
    files = list(raw.glob('*.fits'))
    if len(files) == 0:
        return None

    # search for science files
    modes = []
    for file in files:
        hdr = fits.getheader(file)

        dpr_catg = hdr.get('HIERARCH ESO DPR CATG')
        mode     = hdr.get('HIERARCH ESO INS1 MODE')

        if dpr_catg == 'SCIENCE':
            modes.append(mode)

    modes = np.array(modes)

    n_imaging = (modes == 'NIROBS').sum() + (modes == 'EXT').sum() + \
        (modes == 'DBI').sum() + (modes == 'CI').sum()
    n_pola    = (modes == 'DPI').sum()
    n_spectro = (modes == 'LSS').sum()

    if (n_imaging >= n_pola) and (n_imaging >= n_spectro):
        return 'imaging'
    elif (n_pola >= n_imaging) and (n_pola >= n_spectro):
        return 'pola'
    else:
        return 'spectro'
                

class Dataset:
    '''
    A SPHERE dataset for a given object
    '''

    ######################
    # Constructor
    ######################
    def __init__(self, path, log_level='info'):
        '''
        Initialization code for a SPHERE dataset

        Parameters
        ----------
        path : str
            Path to the SPHERE data
        '''

        if not isinstance(path, str):
            raise ValueError('path must be a string')
        
        # path
        path = Path(path).expanduser().resolve()
        self._path = path

        # configure logging
        logger = logging.getLogger(str(path))
        logger.setLevel(log_level.upper())
        if logger.hasHandlers():
            for hdlr in logger.handlers:
                logger.removeHandler(hdlr)
        
        handler = logging.FileHandler(self._path / 'dataset.log', mode='w', encoding='utf-8')
        formatter = logging.Formatter('%(asctime)s\t%(levelname)8s\t%(message)s')
        formatter.default_msec_format = '%s.%03d'        
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        self._log_level = log_level
        self._handler   = handler
        self._logger    = logger
        
        self._logger.info('Looking for SPHERE data sets at path {}'.format(path))
        
        # list of reductions
        self._IFS_reductions   = []
        self._IRDIS_reductions = []

        # search for data with calibrations downloaded from ESO archive
        xml_files = list(path.glob('*.xml'))
        if len(xml_files) != 0:
            sort_files_from_xml(path, logger=self._logger)

        # directly search for data
        fits_files = list(path.glob('*.fits'))
        if len(fits_files) != 0:
            sort_files_from_fits(path, logger=self._logger)
        
        # recursively look for valid reduction
        self._create_reductions()
        
    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return '<SPHERE datasets: {0} IFS, {1} IRDIS>'.format(len(self._IFS_reductions), len(self._IRDIS_reductions))
    
    ##################################################
    # Properties
    ##################################################
    
    @property
    def reductions(self):
        return self._reductions

    @property
    def IRDIS_reductions(self):
        return self._IRDIS_reductions

    @property
    def IFS_reductions(self):
        return self._IFS_reductions

    ##################################################
    # Generic class methods
    ##################################################
    
    def init_reduction(self):
        '''
        Sort files and frames, perform sanity check
        '''

        for r in self._reductions:
            self._logger.info('Init: {}'.format(str(r)))
            
            r.init_reduction()


    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        for r in self._reductions:
            self._logger.info('Static calibrations: {}'.format(str(r)))
            
            r.create_static_calibrations()

            
    def preprocess_science(self):
        '''
        Clean and collapse images
        '''
        
        for r in self._reductions:
            self._logger.info('Science pre-processing: {}'.format(str(r)))
            
            r.preprocess_science()

            
    def process_science(self):
        '''
        Perform star center, combine cubes into final (x,y,time,lambda)
        cubes, correct anamorphism and scale the images
        '''
        
        for r in self._reductions:
            self._logger.info('Science processing: {}'.format(str(r)))

            r.process_science()

    
    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''
        
        for r in self._reductions:
            print(r)
            self._logger.info('Clean-up: {}'.format(str(r)))

            r.clean()
        
        
    def full_reduction(self):
        '''
        Performs a full reduction of a data set, from the static
        calibrations to the final (x,y,time,lambda) cubes
        '''
        
        for r in self._reductions:
            self._logger.info('###########################################################################')
            self._logger.info('# Full reduction: {}'.format(str(r)))
            self._logger.info('###########################################################################')
            
            r.full_reduction()

    ##################################################
    # Class methods
    ##################################################
        
    def _create_reductions(self):
        '''
        Detect and create valid reductions in path
        '''

        self._logger.info('Create reductions from sorted data')

        wpath = os.walk(self._path)
        for w in wpath:
            subs = w[1]
            if 'raw' in subs:                
                # if directory has a raw/ sub-directory, make sure it
                # has FITS files and that they are from a valid
                # sub-system
                reduction_path = Path(w[0])
                fits_files = list((reduction_path / 'raw').glob('*.fits'))
                if len(fits_files) != 0:
                    hdr = fits.getheader(fits_files[0])
                    try:
                        arm = hdr['HIERARCH ESO SEQ ARM']
                    except KeyError:
                        self._logger.error('No \'HIERARCH ESO SEQ ARM\' keyword in {}'.format(fits_files[0]))
                        continue
                        
                    if arm == 'IRDIS':
                        mode = classify_irdis_dataset(reduction_path, logger=self._logger)

                        # an error occured in dataset classification
                        if mode is None:
                            continue
                        
                        if mode == 'imaging':
                            self._logger.info(' * IRDIS imaging reduction at path {}'.format(reduction_path))
                            reduction  = IRDIS.ImagingReduction(reduction_path, log_level=self._log_level,
                                                                sphere_handler=self._handler)
                        elif mode == 'polar':
                            self._logger.warning('IRDIS DPI not supported yet')
                        elif mode == 'spectro':
                            self._logger.info(' * IRDIS spectro <reduction at path {}'.format(reduction_path))
                            reduction  = IRDIS.SpectroReduction(reduction_path, log_level=self._log_level,
                                                                sphere_handler=self._handler)

                        # save if reduction was successfully created
                        if reduction is not None:
                            self._IRDIS_reductions.append(reduction)
                    elif arm == 'IFS':
                        self._logger.info(' * IFS reduction at path {}'.format(reduction_path))
                        reduction  = IFS.Reduction(reduction_path, log_level=self._log_level, 
                                                   sphere_handler=self._handler)

                        # save if reduction was successfully created
                        if reduction is not None:
                            self._IFS_reductions.append(reduction)
                    else:
                        self._logger.error('Unknown arm {0}'.format(arm))
                        continue

                    # self._logger.info(reduction_path)
                    self._logger.info('   ==> {} files'.format(len(fits_files)))

        # merge all reductions into a single list
        self._reductions = self._IFS_reductions + self._IRDIS_reductions
        
    
