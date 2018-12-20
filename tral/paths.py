# (C) 2014-2015 Elke Schaper

"""
    :synopsis: The system paths module.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import os

# Where are the executables stored? E.g. the tandem repeat detection
# algorithms?
EXEC_DIR = "/usr/local/bin"

# What is the path to the package? E.g. path/to/Tral/tral
PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# The default location of the local config files
CONFIG_DIR = os.path.join(os.path.expanduser('~'), ".tral")

# Where are the data-files, e.g. for pvalue calculation stored?
DATA_DIR = os.path.join(CONFIG_DIR, "data")


def config_file(name):
    ''' Returns the complete path to the config file `name`.

        Returns the complete path to the config file `name`. Searches first in
        `CONFIGDIR`. Otherwise, returns the distribution config file.
    '''

    if os.path.isfile(os.path.join(CONFIG_DIR, name)):
        return os.path.join(CONFIG_DIR, name)
    else:
        return os.path.join(PACKAGE_DIRECTORY, 'data', name)
