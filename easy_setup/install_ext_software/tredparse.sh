#!/bin/bash

# INSTALLING TREDPARSE #####
# TREDPARSE: HLI Short Tandem Repeat (STR) caller 
# https://github.com/humanlongevity/tredparse
# only supported for python > 2.7, python 3 not yet supported
# if tral is saved within a root path one may use to run this script


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


# install virtual environment called "tral2" with python2
virtualenv $TRAL_ENV/python2 -p python2

cd $TRAL_ENV/python2
. bin/activate


## -- TODO possible to use apt within virtenv?

apt-get install python-dev # install header files and static libraries for python2 dev
apt install ipython

######################
### Installation tredparse

pip install tredparse

# Wrapper, sourcing tral2, call tred.py with all arguments "$@"
# with the command "tred" tred.py can be started from everywhere

activateEnv=$TRAL_ENV/python2/bin/activate
echo -e '#!/bin/sh
# wrapper file to easily start tredparse

.' $activateEnv '
tred.py "$@"' > $TRAL_EXT_SOFTWARE/tred


chmod +x $TRAL_EXT_SOFTWARE/tred
cp $TRAL_EXT_SOFTWARE/tred /usr/local/bin # copy wrapper file to system path 
deactivate
cd $HOME

# ######################
# ### Uninstall tredparse

# # symply delete the virtenv in which tred is installed + wrapper file

# rm -rf $TRAL_ENV/python2
# rm $TRAL_EXT_SOFTWARE/tred
# rm /usr/local/bin/tred
