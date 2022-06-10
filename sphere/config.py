import logging
import json

from collections import UserDict

_log = logging.getLogger(__name__)


class Configuration(UserDict):

    ##################################################
    # Constructor
    ##################################################

    def __init__(self, path, logger, config):
        self._path   = path
        self._logger = logger

        super().__init__(config)

    ##################################################
    # dictionary-related functions
    ##################################################

    def __setitem__(self, key, item):
        super().__setitem__(key, item)

        self._logger.debug(f'Saving value {item} for key {key}')

        with open(self._path.root / 'reduction_config.json', 'w') as file:
            file.write(json.dumps(self.data, indent=4))

    def __delitem__(self, key):
        self._logger.error('Configuration keys cannot be modified')

    ##################################################
    # Representation
    ##################################################

    def full_description(self, pad=0):
        repr = ''
        padding = pad*' '

        # parameters
        repr += f'\n{padding}{"Parameter":<30s}Value\n'
        repr += padding + '-'*35 + '\n'
        catgs = ['misc', 'cal', 'preproc', 'center', 'combine', 'clean']
        for catg in catgs:
            keys = [key for key in self if key.startswith(catg)]
            for key in keys:
                repr += f'{padding}{key:<30s}{self[key]}\n'
            repr += padding + '-'*35 + '\n'

        return repr

    def __repr__(self):
        return f'{type(self).__name__}({self.full_description(pad=4)})'

    def __str__(self):
        return self.full_description()
    
    ##################################################
    # Other methods
    ##################################################

    def save(self):
        '''
        Save configuration to reduction directory
        '''
        
        self._logger.info('Saving full config to disk')

        with open(self._path.root / 'reduction_config.json', 'w') as file:
            file.write(json.dumps(self, indent=4))

    def load(self):
        '''
        Load configuration from reduction directory
        '''

        cfg = None
        with open(self._path.root / 'reduction_config.json', 'r') as file:
            cfg = json.load(file)
