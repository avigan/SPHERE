'''
VLT/SPHERE primary module
'''

import os
import glob
import shutil
import math
import pandas as pd
import xml.etree.ElementTree as etree

import vltpf.IRDIS as IRDIS
import vltpf.IFS as IFS

from astropy.io import fits
from astropy.time import Time


def process_mainFiles(mainFiles, files, silent=True):
    '''
    Process top-level file association XML from the ESO archive

    Parameters
    ----------
    mainFiles : etree element
        Main file element

    files : list
        List where the files will be appended

    silent : bool
        Display the activity. Default is True
    '''
    # append file names to the list
    for file in mainFiles:
        fname = file.attrib['name']
        files.append(fname)

        if not silent:
            print(' ==> {0}'.format(fname))


def process_association(tree, files, silent=True):
    '''
    Process file association XML from the ESO archive

    Parameters
    ----------
    tree : etree element
        Associated file element

    files : list
        List where the files will be appended

    silent : bool
        Display the activity. Default is True
    '''
    catg = tree.attrib['category']

    if not silent:
        print(catg)

    # skip unused calibrations
    if (catg == 'IFS_STD_ASTROM') or (catg == 'IFS_STD_PHOT') or \
       (catg == 'IFS_DIST') or (catg == 'IRD_CLI_PHOT') or \
       (catg == 'IRD_DIST'):
        if not silent:
            print(' ==> skipping')
        return

    # process differently mainFiles from associatedFiles
    for elt in tree:
        if elt.tag == 'mainFiles':
            if not silent:
                print('mainFiles')
            process_mainFiles(elt, files)
        elif elt.tag == 'associatedFiles':
            if not silent:
                print('associatedFiles')
            for nelt in elt:
                process_association(nelt, files, silent=silent)


def sort_files_from_xml(path, silent=True):
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
    silent : bool
        Display some status of the execution. Default is to be silent

    '''

    xml_files = glob.glob(path+'*.xml')

    print('Sort data based on XML files (ESO automated calibration selection)')
    print(' ==> {0} XML files\n'.format(len(xml_files)))
    
    # sort files
    for file in xml_files:
        tree = etree.parse(file)
        root = tree.getroot()            

        print(os.path.basename(file))
        
        # process only IFS and IRDIS science data
        catg = root.attrib['category']
        if (catg.find('ACQUISITION') != -1):
            continue

        # get target name from first mainFile element
        scifiles = root.find('mainFiles')
        filename = scifiles[0].attrib['name']
        
        # Mac OS X replaces : by _ in file names...
        if not os.path.exists(path+filename+'.fits'):
            filename = filename.replace(':', '_')
        
        if not os.path.exists(path+filename+'.fits'):
            print(' ==> file {} does not exsist. Skipping'.format(filename))
            continue

        hdr = fits.getheader(path+filename+'.fits')
        
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
                if arm == 'IRDIS':
                    instrument = 'IRDIS'
                elif arm == 'IFS':
                    instrument = 'IFS'                
                else:
                    raise NameError('Unknown arm {0}'.format(arm))
            except NameError:
                continue

        # get files
        files = []
        process_association(root, files, silent=True)

        # target path
        directory = '{0}_id={1}'.format(target, obs_id)
        directory = '_'.join(directory.split())
        target_path = os.path.join(path, directory, night, instrument, 'raw')
        if not os.path.exists(target_path):
            os.makedirs(target_path)

        # copy files
        for f in files:
            fpath  = os.path.join(path, f+'.fits')
            
            # Mac OS X replaces : by _ in file names...
            if not os.path.exists(fpath):
                tmp = f.replace(':', '_')
                fpath  = os.path.join(path, tmp+'.fits')
            
            # check if file actually exists
            if not os.path.exists(fpath):
                print(' ==> file {} does not exist. Skipping.'.format(fpath))
                continue
            
            # copy if needed
            nfpath = os.path.join(target_path, f+'.fits')
            if not os.path.exists(nfpath):
                shutil.copy(fpath, nfpath)

        # print status
        if not silent:
            print('{0} - id={1}'.format(target, obs_id))
            print(' ==> found {0} files'.format(len(files)))
            print(' ==> copied to {0}'.format(target_path))
            print()

    # move all files
    path_new = os.path.join(path, 'all_files')
    if not os.path.exists(path_new):
        os.makedirs(path_new)

    files = []
    files.extend(glob.glob(os.path.join(path+'*.fits')))
    files.extend(glob.glob(os.path.join(path+'*.xml')))
    files.extend(glob.glob(os.path.join(path+'*.txt')))        

    if len(files) != 0:
        for file in files:
            shutil.move(file, path_new)


def sort_files_from_fits(path, silent=True):
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
    silent : bool
        Display some status of the execution. Default is to be silent

    '''

    fits_files = glob.glob(path+'*.fits')

    print('Sort data based on FITS files')
    print(' ==> {0} FITS files\n'.format(len(fits_files)))

    # sort files
    for file in fits_files:
        # target and arm
        hdr = fits.getheader(file)

        try:
            target = hdr['HIERARCH ESO OBS NAME']
            obs_id = hdr['HIERARCH ESO OBS ID']
            dpr_type = hdr['HIERARCH ESO DPR TYPE']
        except:
            continue

        if dpr_type == 'OBJECT,AO':
            instrument = 'SPARTA'
        else:
            try:
                arm = hdr['HIERARCH ESO SEQ ARM']
                if arm == 'IRDIS':
                    instrument = 'IRDIS'
                elif arm == 'IFS':
                    instrument = 'IFS'                
                else:
                    raise NameError('Unknown arm {0}'.format(arm))
            except:
                continue

        # target path
        target_path = os.path.join(path, instrument, 'raw')
        if not os.path.exists(target_path):
            os.makedirs(target_path)

        # move file
        nfile = os.path.join(target_path, os.path.basename(file))
        shutil.move(file, nfile)

        # print status
        if not silent:
            print('{0} - id={1}'.format(target, obs_id))
            print(' ==> copied to {0}'.format(target_path))
            print()

    # move all files
    path_new = os.path.join(path, 'unsorted_files')
    if not os.path.exists(path_new):
        os.makedirs(path_new)

    files = []
    files.extend(glob.glob(os.path.join(path+'*.fits')))
    files.extend(glob.glob(os.path.join(path+'*.txt')))        

    if len(files) != 0:
        for file in files:
            shutil.move(file, path_new)


class Dataset:
    '''
    A SPHERE dataset for a given object
    '''

    ######################
    # Constructor
    ######################
    def __init__(self, path):
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
        path = os.path.expanduser(os.path.join(path, ''))
        self._path = path

        # list of reductions
        self._IFS_reductions   = []
        self._IRDIS_reductions = []

        # search for data with calibrations downloaded from ESO archive
        xml_files = glob.glob(os.path.join(path, '*.xml'))
        if len(xml_files) != 0:
            sort_files_from_xml(path)

        # directly search for data
        fits_files = glob.glob(os.path.join(path, '*.fits'))
        if len(fits_files) != 0:
            sort_files_from_fits(path)
        
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
            print()
            print('*')
            print('* Initialization of {0} reduction at path {1}'.format(r.instrument, r.path))
            print('*')
            print()
            
            r.init_reduction()


    def create_static_calibrations(self):
        '''
        Create static calibrations with esorex
        '''

        for r in self._reductions:
            print()
            print('*')
            print('* Static calibrations for {0} at path {1}'.format(r.instrument, r.path))
            print('*')
            print()
            
            r.create_static_calibrations()

            
    def preprocess_science(self):
        '''
        Clean and collapse images
        '''
        
        for r in self._reductions:
            print()
            print('*')
            print('* Pre-process data for {0} at path {1}'.format(r.instrument, r.path))
            print('*')
            print()
            
            r.preprocess_science()

            
    def process_science(self):
        '''
        Perform star center, combine cubes into final (x,y,time,lambda)
        cubes, correct anamorphism and scale the images
        '''
        
        for r in self._reductions:
            print()
            print('*')
            print('* Process data for {0} at path {1}'.format(r.instrument, r.path))
            print('*')
            print()

            r.process_science()

    
    def clean(self):
        '''
        Clean the reduction directory, leaving only the raw and products
        sub-directory
        '''
        
        for r in self._reductions:
            print()
            print('*')
            print('* Clean {0} reduction at path {1}'.format(r.instrument, r.path))
            print('*')
            print()

            r.clean()
        
        
    def full_reduction(self):
        '''
        Performs a full reduction of a data set, from the static
        calibrations to the final (x,y,time,lambda) cubes
        '''
        
        for r in self._reductions:
            print()
            print('*')
            print('* Full {0} reduction at path {1}'.format(r.instrument, r.path))
            print('*')
            print()
            
            r.full_reduction()

            
    ##################################################
    # Class methods
    ##################################################
        
    def _create_reductions(self):
        '''
        Detect and create valid reductions in path
        '''

        print('Create reductions from available data')

        wpath = os.walk(self._path)
        for w in wpath:
            subs = w[1]
            if 'raw' in subs:                
                # if directory has a raw/ sub-directory, make sure it
                # has FITS files and that they are from a valid
                # sub-system
                reduction_path = w[0]                
                fits_files = glob.glob(os.path.join(reduction_path, 'raw', '*.fits'))
                if len(fits_files) != 0:
                    hdr = fits.getheader(fits_files[0])
                    try:
                        arm = hdr['HIERARCH ESO SEQ ARM']
                        if arm == 'IRDIS':
                            instrument = 'IRDIS'
                            reduction  = IRDIS.ImagingReduction(reduction_path)
                            self._IRDIS_reductions.append(reduction)
                        elif arm == 'IFS':
                            instrument = 'IFS'
                            reduction  = IFS.Reduction(reduction_path)
                            self._IFS_reductions.append(reduction)
                        else:
                            raise NameError('Unknown arm {0}'.format(arm))
                    except:
                        continue

                    print(reduction_path)
                    print('  ==> {0}, {1} files'.format(instrument, len(fits_files)))
                    print()

        # merge all reductions into a single list
        self._reductions = self._IFS_reductions + self._IRDIS_reductions
        
    
