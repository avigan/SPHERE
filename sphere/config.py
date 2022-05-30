import logging
import json

_log = logging.getLogger(__name__)

class Configuration(dict):

    ##################################################
    # Constructor
    ##################################################

    def __init__(self, instrument, path, logger, inpt={}):
        super().__init__(inpt)

        self._path   = path
        self._logger = logger

    def __setitem__(self, key, item):
        # self.__dict__[key] = item
        super().__setitem__(key, item)

        self._logger.debug(f'Saving value {item} for key {key}')

        with open(self._path.root / 'reduction_config.json', 'w') as file:
            file.write(json.dumps(self, indent=4))
    
    ##################################################
    # Representation
    ##################################################

    def __repr__(self):
        repr = ''

        # parameters
        repr += f'\n{"Parameter":<30s}Value\n'
        repr += '-'*35 + '\n'
        catgs = ['misc', 'cal', 'preproc', 'center', 'combine', 'clean']
        for catg in catgs:
            keys  = [key for key in self if key.startswith(catg)]
            for key in keys:
                repr += f'{key:<30s}{self[key]}\n'
            repr += '-'*35 + '\n'
        
        return repr

    def __format__(self):
        return self.__repr__()
