__author__    = 'Arthur Vigan'
__copyright__ = 'Copyright (C) 2017-2020 Arthur Vigan'
__license__   = 'MIT'
__version__   = '1.2'

import logging
import enum

#
# recipe execution status
#
class RecipeStatus(enum.IntEnum):
    NOTSET     = -1
    SUCCESS    =  0
    ERROR      =  1
    FATAL      =  2

NOTSET     = RecipeStatus.NOTSET
SUCCESS    = RecipeStatus.SUCCESS
ERROR      = RecipeStatus.ERROR
FATAL      = RecipeStatus.FATAL

#
# reduction execution status
#
class ReductionStatus(enum.IntEnum):
    COMPLETE   =   0
    INIT       = -10
    INCOMPLETE = -20

COMPLETE   = ReductionStatus.COMPLETE
INIT       = ReductionStatus.INIT
INCOMPLETE = ReductionStatus.INCOMPLETE

#
# logging
#

# define logging format for module
logging.basicConfig(format='[%(levelname)8s] %(message)s')
_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)
_log.info('vlt-sphere init')

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
