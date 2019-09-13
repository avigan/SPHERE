__author__ = 'avigan'
__copyright__ = 'Copyright (C) 2017-2019 Arthur Vigan'
__license__ = 'MIT'

import logging

# define logging format for module
logging.basicConfig(format='[%(levelname)8s] %(message)s')
_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)
_log.info('VLTPF init')

# reduction status
SUCCESS    =  0
ERROR      =  1
INIT       = -1
INCOMPLETE = -2


# log level
def set_loglevel(level):
    '''
    Set the logging level for the module

    Parameters
    ----------
    level : {"notset", "debug", "info", "warning", "error", "critical"}
        The log level of the handler
    '''
    
    _log.setLevel(level.upper())
