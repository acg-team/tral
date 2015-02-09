import os

# Where are the executables stored? E.g. the tandem repeat detection algorithms?
EXEC_DIR = "/usr/local/bin"

# What is the path to the package? E.g. path/to/Tral/tral
PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# The default location of the local config files
CONFIG_DIR = os.path.join(os.path.expanduser('~'), ".tral")

# Where are the data-files, e.g. for pValue calculation stored?
DATA_DIR = os.path.join(CONFIG_DIR, "data")

def config():

    if os.path.isfile(CONFIG_DIR, "config.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "config.ini")
    else:
        return os.path.join(DATAROOT, "config.ini")

def config_spec():

    if os.path.isfile(CONFIG_DIR, "spec.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "spec.ini")
    else:
        return os.path.join(DATAROOT, "spec.ini")

def logging_spec():

    if os.path.isfile(CONFIG_DIR, "logging.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "logging.ini")
    else:
        return os.path.join(DATAROOT, "logging.ini")
