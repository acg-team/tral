#!/bin/bash

# INSTALLING MAFFT #####
# MAFFT: Multiple alignment program for amino acid or nucleotide sequences

# MAFFT is a multiple sequence alignment program for unix-like operating systems.
# It offers a range of multiple alignment methods, L-INS-i (accurate; for alignment of <∼200 sequences), FFT-NS-2 (fast; for alignment of <∼30,000 sequences), etc. 
# http://hmmer.org/documentation.html


# The mafft-6.6xx-with-extensions-src.tgz package has the 'extensions' directory, which consists of the codes from:
# (1) the Vienna RNA package 
# (2) MXSCARNA
# (3) ProbConsRNA
# These are distributed under different licenses: https://mafft.cbrc.jp/alignment/software/license66.txt



######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Download and Installation MAFFT

latestVer=$(wget -qO- https://mafft.cbrc.jp/alignment/software/linux.html |
    egrep amd64.deb |                  # only grep deb version
    sed -n 's/.*href="\([^"]*\).*/\1/p')

dpkg -s mafft 2>/dev/null >/dev/null && echo "MAFFT is already installed " || # test if mafft already installed
    {
    # if not installed get newest version for linux system
    if [ ! -f $TRAL_EXT_SOFTWARE/$latestVer ]; then # test if latest version already downloaded
        wget https://mafft.cbrc.jp/alignment/software/$latestVer -P $TRAL_EXT_SOFTWARE # download
    fi
    dpkg -i $TRAL_EXT_SOFTWARE/$latestVer # install latest version
    }

# mafft is now accessible with the command "mafft"

######################
### Uninstallation of MAFFT (default paths!)

# dpkg --remove mafft
# rm -rf $TRAL_EXT_SOFTWARE/mafft*
