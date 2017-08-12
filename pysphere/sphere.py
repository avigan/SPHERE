'''
VLT/SPHERE primary module
'''

import os
import glob
import shutil
import numpy as np
import pandas as pd

from astropy.io import fits


class IRDISData:
    pass


class IFSData:
    pass


class SPHEREDataset:
    '''
    A SPHERE dataset for a given object
    '''

    rootpath = None

    # IFS
    ifs_data   = False
    ifs_sorted = False
    
    # IRD
    ird_data   = False
    ird_sorted = False
    
    ######################
    # Constructor
    ######################
    def __init__(self, rootpath):
        '''
        Initialization code for a SPHERE dataset

        Parameters
        ----------
        rootpath : str
            Root path to the SPHERE data set
        '''

        if not isinstance(rootpath, str):
            raise ValueError('rootpath must be a string')
        
        self.rootpath = os.path.expanduser(rootpath)

        # check if data is already sorted or not
        ipath = os.path.join(rootpath, 'IFS', 'raw')
        if len(glob.glob(ipath+'/*.fits')) != 0:
            self.ifs_data = True
            
        ipath = os.path.join(rootpath, 'IRD', 'raw')
        if len(glob.glob(ipath+'/*.fits')) != 0:
            self.ird_data = True

    
    def sort_files(self):
        '''
        Sort the raw files in the rootpath directory        
        '''

        files = glob.glob(self.rootpath+'*.fits')

        # check that we have some files
        if len(files) == 0:
            raise ValueError('No raw FITS files in rootpath directory')

        print('Found {0} FITS files in {1}'.format(len(files), self.rootpath))
        
        # sort them by sub-system
        for f in files:
            hdu = fits.open(f)
            subsystem = hdu[0].header['HIERARCH ESO SEQ ARM']
            hdu.close()

            # define instrument short name
            if subsystem == 'IFS':
                inst = 'IFS'
            elif subsystem == 'IRDIS':
                inst = 'IRD'
            elif subsystem == 'ZIMPOL':
                inst = 'ZIM'

            # move file
            ipath = os.path.join(self.rootpath, inst, 'raw')
            if not os.path.exists(ipath):
                os.makedirs(ipath)
            shutil.move(f, ipath)
                    





root = '/Users/avigan/data/pySPHERE-test/'

sphere_ds = SPHEREDataset(root)
# sphere_ds.sort_files()
