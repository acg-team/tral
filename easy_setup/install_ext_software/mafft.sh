#!/bin/bash

# INSTALLING MAFFT #####
# MAFFT: Multiple alignment program for amino acid or nucleotide sequences
# http://hmmer.org/documentation.html
# 


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Download and Installion MAFFT

# get newest version for linux system
latestVer=$(wget -qO- https://mafft.cbrc.jp/alignment/software/linux.html |
            egrep amd64.deb |                  # only grep deb version
            sed -n 's/.*href="\([^"]*\).*/\1/p')

sudo wget https://mafft.cbrc.jp/alignment/software/$latestVer -P $TRAL_EXT_SOFTWARE # download
sudo dpkg -i $TRAL_EXT_SOFTWARE/$latestVer


######################
### Uninstall MAFFT (default paths!)
# TODO -- check if this uninstalls also if package installed within the $TRAL_EXT_SOFTWARE

# sudo dpkg --remove mafft
# sudo rm -rf $TRAL_EXT_SOFTWARE/mafft*
