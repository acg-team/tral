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

latestVer=$(wget -qO- https://mafft.cbrc.jp/alignment/software/linux.html |
    egrep amd64.deb |                  # only grep deb version
    sed -n 's/.*href="\([^"]*\).*/\1/p')

dpkg -s mafft 2>/dev/null >/dev/null && echo "MAFFT is already installed " || # test if mafft already installed
    {
    # if not installed get newest version for linux system
    wget https://mafft.cbrc.jp/alignment/software/$latestVer -P $TRAL_EXT_SOFTWARE # download
    dpkg -i $latestVer
    }



######################
### Uninstall MAFFT (default paths!)

# dpkg --remove mafft
# rm -rf $TRAL_EXT_SOFTWARE/mafft*
