'''
VLT/SPHERE primary module
'''

import os
import glob
import shutil
import pandas as pd
import xml.etree.ElementTree as etree

import pysphere.IRDIS as IRDIS
import pysphere.IFS as IFS

from astropy.io import fits


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
        self._path = os.path.expanduser(path)

        # list of reductions
        self._IFS_reductions   = []
        self._IRDIS_reductions = []

        # search for data with calibrations downloaded from ESO archive
        xml_files = glob.glob(path+'*.xml')
        if len(xml_files) != 0:
            print('Searching for for data with calibrations downloaded from ESO archive')
            self.sort_files_from_archive()

        # recursively look for valid reduction
        wpath = os.walk(path)
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
    
    def sort_files_from_archive(self, silent=True):
        '''Sort files downloaded from the ESO archive with associated raw
           calibrations

        When choosing this option from the archive, the data package
        includes xml files that provide the association between the
        science files and the calibration files. The method uses these
        xml files to copy the science and associated raw files into
        dedicated directories.

        For a given executed OB, the method will extract the target
        name from the OBS.NAME keyword and the OB ID from the OBS.ID
        keyword. It will create a directory {OBD.NAME}_id={OBS.ID} and
        copy the FITS files. Note that any space in the name will be
        replaced by an underscore to avoid subsequent issues with
        esorex in the data reduction.

        Parameters
        ----------
        silent : bool
            Display some status of the execution. Default is to be silent
        '''
        path = self._path

        xml_files = glob.glob(path+'*.xml')

        if len(xml_files) == 0:
            print('This path does not appear to contain a dataset downloaded from ' + 
                  'the ESO archive with associated calibrations. Skipping.')
            return

        # sort files
        for file in xml_files:
            tree = etree.parse(file)
            root = tree.getroot()            

            # process only IFS and IRDIS science data
            catg = root.attrib['category']
            if (catg.find('ACQUISITION') != -1):
                continue
            
            # get target name from first mainFile element
            scifiles = root.find('mainFiles')
            filename = scifiles[0].attrib['name']

            # target and arm
            hdu = fits.open(path+filename+'.fits')

            target = hdu[0].header['HIERARCH ESO OBS NAME']
            obs_id = hdu[0].header['HIERARCH ESO OBS ID']
            if catg == 'SCIENCE_OBJECT_AO':
                instrument = 'SPARTA'
            else:
                try:
                    arm = hdu[0].header['HIERARCH ESO SEQ ARM']
                    if arm == 'IRDIS':
                        instrument = 'IRDIS'
                    elif arm == 'IFS':
                        instrument = 'IFS'                
                    else:
                        raise NameError('Unknown arm {0}'.format(arm))
                except:
                    continue

            # get files
            files = []
            process_association(root, files, silent=True)

            # target path
            directory = '{0}_id={1}'.format(target, obs_id)
            directory = '_'.join(directory.split())
            target_path = os.path.join(path, directory, instrument, 'raw')
            if not os.path.exists(target_path):
                os.makedirs(target_path)

            # copy files
            for f in files:
                file = os.path.join(path, f+'.fits')
                nfile = os.path.join(target_path, f+'.fits')

                # copy if needed
                if not os.path.exists(nfile):
                    shutil.copy(file, nfile)

            # print status
            if not silent:
                print('{0} - id={1}'.format(target, obs_id))
                print(' ==> found {0} files'.format(len(files)))
                print(' ==> copied to {0}'.format(target_path))
                print()

        # move all files
        path_new = os.path.join(path, 'files')
        if not os.path.exists(path_new):
            os.makedirs(path_new)
            
        files = []
        files.extend(glob.glob(os.path.join(path+'*.fits')))
        files.extend(glob.glob(os.path.join(path+'*.xml')))
        files.extend(glob.glob(os.path.join(path+'*.txt')))        

        if len(files) != 0:
            for file in files:
                shutil.move(file, path_new)
