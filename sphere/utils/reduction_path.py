import shutil
import logging

from pathlib import Path

_log = logging.getLogger(__name__)


class ReductionPath(object):
    '''
    Reduction path class

    Provides easy access to all the path components of a reduction.
    '''

    ##################################################
    # Constructor
    ##################################################

    def __init__(self, path):
        self._root = Path(path).expanduser()

        # update all subpaths
        self._raw      = self._root / 'raw'
        self._calib    = self._root / 'calib'
        self._sof      = self._root / 'sof'
        self._tmp      = self._root / 'tmp'
        self._preproc  = self._root / 'preproc'
        self._products = self._root / 'products'

    ##################################################
    # Representation
    ##################################################

    def __repr__(self):
        return str(self._root)

    ##################################################
    # Properties
    ##################################################

    @property
    def root(self):
        return self._root

    @root.setter
    def root(self, path):
        self._root = Path(path).expanduser()

        # update all subpaths
        self._raw      = self._root / 'raw'
        self._calib    = self._root / 'calib'
        self._sof      = self._root / 'sof'
        self._tmp      = self._root / 'tmp'
        self._preproc  = self._root / 'preproc'
        self._products = self._root / 'products'

    @property
    def raw(self):
        # create sub-directory if needed
        if not self._raw.exists():
            self._raw.mkdir(exist_ok=True)

        return self._raw

    @property
    def calib(self):
        # create sub-directory if needed
        if not self._calib.exists():
            self._calib.mkdir(exist_ok=True)

        return self._calib

    @property
    def sof(self):
        # create sub-directory if needed
        if not self._sof.exists():
            self._sof.mkdir(exist_ok=True)

        return self._sof

    @property
    def tmp(self):
        # create sub-directory if needed
        if not self._tmp.exists():
            self._tmp.mkdir(exist_ok=True)

        return self._tmp

    @property
    def preproc(self):
        # create sub-directory if needed
        if not self._preproc.exists():
            self._preproc.mkdir(exist_ok=True)

        return self._preproc

    @property
    def products(self):
        # create sub-directory if needed
        if not self._products.exists():
            self._products.mkdir(exist_ok=True)

        return self._products

    ##################################################
    # Properties
    ##################################################

    def remove(self, delete_raw=False, delete_products=False, logger=_log):
        '''
        Remove sub-directories

        Parameters
        ----------
        delete_raw : bool
            Delete raw data. Default is False

        delete_products : bool
            Delete science products. Default is False

        logger : logHandler object
            Log handler for the reduction. Default is root logger
        '''

        # tmp
        if self._tmp.exists():
            logger.debug(f'> remove {self._tmp}')
            shutil.rmtree(self._tmp, ignore_errors=True)

        # sof
        if self._sof.exists():
            logger.debug(f'> remove {self._sof}')
            shutil.rmtree(self._sof, ignore_errors=True)

        # calib
        if self._calib.exists():
            logger.debug(f'> remove {self._calib}')
            shutil.rmtree(self._calib, ignore_errors=True)

        # preproc
        if self._preproc.exists():
            logger.debug(f'> remove {self._preproc}')
            shutil.rmtree(self._preproc, ignore_errors=True)

        # raw
        if delete_raw:
            if self._raw.exists():
                logger.debug(f'> remove {self._raw}')
                logger.warning('   ==> delete raw files')
                shutil.rmtree(self._raw, ignore_errors=True)

        # products
        if delete_products:
            if self._products.exists():
                logger.debug(f'> remove {self._products}')
                logger.warning('   ==> delete products')
                shutil.rmtree(self._products, ignore_errors=True)
