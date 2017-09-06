'''
VLT/SPHERE primary module
'''

import os
import glob
import shutil
import pandas as pd
import xml.etree.ElementTree as etree

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
            raise ValueError('rootpath must be a string')

        # path
        self._path = os.path.expanduser(path)

        # list of reductions
        self._IFS_reductions   = []
        self._IRDIS_reductions = []

        # search for available data
        print(os.listdir(path))
        
    # def sort_files(self):
    #     '''
    #     Sort the raw files in the rootpath directory        
    #     '''

    #     files = glob.glob(self.rootpath+'*.fits')

    #     # check that we have some files
    #     if len(files) == 0:
    #         raise ValueError('No raw FITS files in rootpath directory')

    #     print('Found {0} FITS files in {1}'.format(len(files), self.rootpath))
        
    #     # sort them by sub-system
    #     for f in files:
    #         hdu = fits.open(f)
    #         subsystem = hdu[0].header['HIERARCH ESO SEQ ARM']
    #         hdu.close()

    #         # define instrument short name
    #         if subsystem == 'IFS':
    #             inst = 'IFS'
    #         elif subsystem == 'IRDIS':
    #             inst = 'IRD'
    #         elif subsystem == 'ZIMPOL':
    #             inst = 'ZIM'

    #         # move file
    #         ipath = os.path.join(self.rootpath, inst, 'raw')
    #         if not os.path.exists(ipath):
    #             os.makedirs(ipath)
    #         shutil.move(f, ipath)

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

        xfiles = glob.glob(path+'*.xml')

        if len(xfiles) == 0:
            print('This path does not appear to contain a dataset downloaded from ' + 
                  'the ESO archive with associated calibrations. Skipping.')
            return

        for xfile in xfiles:
            tree = etree.parse(xfile)
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
