#!/bin/bash

# INSTALLING ALF #####

######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file

[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE/bin:$PATH"* ]] && export PATH="$TRAL_EXT_SOFTWARE/bin:$PATH"


######################
### Download and Installion ALF

if [ ! -d $TRAL_EXT_SOFTWARE/ALF_standalone ]; then # test if not already in directory
    LINK_ALF=http://abacus.gene.ucl.ac.uk/daniel/alf/ALF_standalone.tar.gz
    wget $LINK_ALF -P $TRAL_EXT_SOFTWARE    # download
    tar -xvzf $TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz -C $TRAL_EXT_SOFTWARE
fi

rm -rf $TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz
(cd $TRAL_EXT_SOFTWARE/ALF_standalone && $TRAL_EXT_SOFTWARE/ALF_standalone/install.sh) # installation of ALF


######################
### Uninstall ALF (default paths!)

# rm -rf /usr/local/bin/darwin*
# rm -rf /usr/local/bin/alfsim
# rm -rf /usr/local/share/alfdarwin


