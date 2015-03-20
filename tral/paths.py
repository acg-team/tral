import os

# Where are the executables stored? E.g. the tandem repeat detection algorithms?
EXEC_DIR = "/usr/local/bin"

# What is the path to the package? E.g. path/to/Tral/tral
PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# The default location of the local config files
CONFIG_DIR = os.path.join(os.path.expanduser('~'), ".tral")

# Where are the data-files, e.g. for pValue calculation stored?
DATA_DIR = os.path.join(CONFIG_DIR, "data")

def config_file(name):

    if os.path.isfile(os.path.join(CONFIG_DIR, name)):
        return os.path.join(CONFIG_DIR, name)
    else:
        return os.path.join(PACKAGE_DIRECTORY, 'data', name)
