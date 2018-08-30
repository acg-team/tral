#!/bin/bash

# INSTALLING ALF #####

# (!) this script will install ALF within the $TRAL_EXT_SOFTWARE directory

######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file

[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE/bin:$PATH"* ]] && export PATH="$TRAL_EXT_SOFTWARE/bin:$PATH"


######################
### Download and Installion ALF

LINK_ALF=http://abacus.gene.ucl.ac.uk/daniel/alf/ALF_standalone.tar.gz
sudo wget $LINK_ALF -P $TRAL_EXT_SOFTWARE    # download
sudo tar -xvzf $TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz -C $TRAL_EXT_SOFTWARE
sudo rm -rf $TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz
(cd $TRAL_EXT_SOFTWARE/ALF_standalone && sudo $TRAL_EXT_SOFTWARE/ALF_standalone/install.sh $TRAL_EXT_SOFTWARE) # installation of ALF


######################
### Uninstall ALF (default paths!)

# rm -rf /usr/local/bin/darwin*
# rm -rf /usr/local/bin/alfsim
# rm -rf /usr/local/share/alfdarwin


