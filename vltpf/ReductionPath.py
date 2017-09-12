import os


class Path(object):
    '''
    Reduction path class

    Provides easy access to all the path components of a reduction.
    '''

    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, path):
        self._root = path

        # update all subpaths
        self._raw = os.path.join(self._root, 'raw/')
        self._calib = os.path.join(self._root, 'calib/')
        self._sof = os.path.join(self._root, 'sof/')
        self._tmp = os.path.join(self._root, 'tmp/')
        self._preproc = os.path.join(self._root, 'preproc/')
        self._products = os.path.join(self._root, 'products/')

        # create directories
        self.create_subdirectories()

    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return self._root
    
    ##################################################
    # Properties
    ##################################################
    
    @property
    def root(self):
        return self._root

    @root.setter
    def root(self, path):
        self._root = os.path.expanduser(path)

        # update all subpaths
        self._raw = os.path.join(self._root, 'raw/')
        self._calib = os.path.join(self._root, 'calib/')
        self._sof = os.path.join(self._root, 'sof/')
        self._tmp = os.path.join(self._root, 'tmp/')
        self._preproc = os.path.join(self._root, 'preproc/')
        self._products = os.path.join(self._root, 'products/')

        # create directories
        self.create_subdirectories()

    @property
    def raw(self):
        return self._raw
    
    @property
    def calib(self):
        return self._calib
    
    @property
    def sof(self):
        return self._sof
    
    @property
    def tmp(self):
        return self._tmp
    
    @property
    def preproc(self):
        return self._preproc
    
    @property
    def products(self):
        return self._products

    ##################################################
    # Methods
    ##################################################
    
    def create_subdirectories(self):
        # create sub-directories if needed
        if not os.path.exists(self._raw):
            os.makedirs(self._raw)

        if not os.path.exists(self._calib):
            os.makedirs(self._calib)

        if not os.path.exists(self._sof):
            os.makedirs(self._sof)

        if not os.path.exists(self._tmp):
            os.makedirs(self._tmp)

        if not os.path.exists(self._preproc):
            os.makedirs(self._preproc)

        if not os.path.exists(self._products):
            os.makedirs(self._products)

