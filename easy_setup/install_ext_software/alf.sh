#!/bin/bash

# INSTALLING ALF #####

######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file

[[ ":$PATH:" != *"$TRAL_SOFTWARE/bin:$PATH"* ]] && export PATH="$TRAL_SOFTWARE/bin:$PATH"


######################
### Download and Installion ALF
LINK_ALF=http://abacus.gene.ucl.ac.uk/daniel/alf/ALF_standalone.tar.gz
sudo wget $LINK_ALF -P $TRAL_EXT
sudo tar -xvzf $TRAL_EXT/ALF_standalone.tar.gz -C $TRAL_SOFTWARE
sudo rm -rf $TRAL_EXT/ALF_standalone.tar.gz
(cd $TRAL_SOFTWARE/ALF_standalone && sudo $TRAL_SOFTWARE/ALF_standalone/install.sh $TRAL_SOFTWARE) # installation of ALF

alfsim