import os

# Where are the executables stored? E.g. the tandem repeat detection algorithms?
EXECROOT = "/usr/local/bin"

# What is the path to the package? E.g. path/to/Tral/tral
PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# Where are the data-files, e.g. for pValue calculation stored?
DATAROOT = os.path.join(PACKAGE_DIRECTORY, "data")

def config():

    if os.path.isfile(os.path.join(os.path.expanduser('~'), ".tral", "config.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "config.ini")
    else:
        return os.path.join(DATAROOT, "config.ini")

def config_spec():

    if os.path.isfile(os.path.join(os.path.expanduser('~'), ".tral", "spec.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "spec.ini")
    else:
        return os.path.join(DATAROOT, "spec.ini")

def configuration():

    if os.path.isfile(os.path.join(os.path.expanduser('~'), ".tral", "logging.ini")):
        return os.path.join(os.path.expanduser('~'), ".tral", "logging.ini")
    else:
        return os.path.join(DATAROOT, "logging.ini")
