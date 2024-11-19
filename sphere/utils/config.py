import logging
import configparser

from collections import UserDict

_log = logging.getLogger(__name__)


class Configuration(UserDict):

    ##################################################
    # Constructor
    ##################################################

    def __init__(self, path, logger, config):
        self._file   = path.root / 'reduction_config.ini'
        self._logger = logger

        # initialize internal dict with user-provided configuration
        self.data = config

        ##################################################
    # dictionary-related functions
    ##################################################

    def __setitem__(self, key, item):
        super().__setitem__(key, item)

        self._logger.debug(f'Saving value {item} for key {key}')

        self.save()

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

        self._logger.debug('Saving full config to disk')

        with open(self._file, 'w') as file:
            file.write('[default]\n\n')

            catgs = ['misc', 'cal', 'preproc', 'center', 'combine', 'clean']
            for catg in catgs:
                keys = [key for key in self if key.startswith(catg)]
                for key in keys:
                    file.write(f'{key:<30s} = {self[key]}\n')
                file.write('\n')

    def load(self):
        '''
        Load configuration from reduction directory
        '''

        if self._file.exists():
            self._logger.info('Load existing configuration file')

            try:
                cfgparser = configparser.ConfigParser()
                cfgparser.read(self._file)
            except configparser.MissingSectionHeaderError:
                # add section if it was missing
                file = open(self._file, 'r')
                lines = file.readlines()

                with open(self._file, 'w') as file:
                    file.write('[default]\n\n')
                    file.writelines(lines)
            finally:
                cfgparser = configparser.ConfigParser()
                cfgparser.read(self._file)

                for section in cfgparser.sections():
                    items = dict(cfgparser.items(section))
                    for key, value in items.items():
                        try:
                            val = eval(value)
                        except NameError:
                            val = value

                        self.data[key] = val
        else:
            self.save()

    def load_from_file(self, filepath):
        '''
        Load configuration from provided file

        Parameters
        ----------
        filepath : str
            Path of the configuration file
        '''

        if filepath.exists():
            self._logger.info(f'Load configuration file at path {filepath}')

            try:
                cfgparser = configparser.ConfigParser()
                cfgparser.read(filepath)
            except configparser.MissingSectionHeaderError:
                # add section if it was missing
                file = open(filepath, 'r')
                lines = file.readlines()

                with open(filepath, 'w') as file:
                    file.write('[default]\n\n')
                    file.writelines(lines)
            finally:
                cfgparser = configparser.ConfigParser()
                cfgparser.read(filepath)

                for section in cfgparser.sections():
                    items = dict(cfgparser.items(section))
                    for key, value in items.items():
                        try:
                            val = eval(value)
                        except NameError:
                            val = value

                        self.data[key] = val

