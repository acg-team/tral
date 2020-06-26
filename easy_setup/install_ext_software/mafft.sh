#!/usr/bin/env bash

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

shopt -s nocasematch # making comparisons case-insensitive
set -euo pipefail # exit on errors and undefined vars

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit 1
}

# ######################
# ### Download and Installation MAFFT
latestVer="mafft-7.464-without-extensions" # CHANGE MAFFT VERSION HERE
mafftVer=${latestVer}-src.tgz

echo "---------------------------------------------------------------------"
echo "Downloading:" $latestVer
echo "Change version in tral/easy_setup/install_ext_software/mafft.sh to"
echo "the latest version from https://mafft.cbrc.jp/alignment/software/source.html"
echo "---------------------------------------------------------------------"

sleep 5

if [ ! -d "$TRAL_EXT_SOFTWARE/$mafftVer" ]; then
    {
        wget "https://mafft.cbrc.jp/alignment/software/"$mafftVer -P "$TRAL_EXT_SOFTWARE" || {  # download
            echo "Was not able to download mafft."
            exit 1
        }
        tar -xvzf "$TRAL_EXT_SOFTWARE/$mafftVer" -C "$TRAL_EXT_SOFTWARE"
        sed -i "s#PREFIX = /usr/local#PREFIX = \"$INSTALLATION_PATH\"#" "$TRAL_EXT_SOFTWARE/$latestVer/core/Makefile" # change default installation path in Makefile
        sed -i "s#BINDIR = \$(PREFIX)/bin#BINDIR = \$(PREFIX)#" "$TRAL_EXT_SOFTWARE/$latestVer/core/Makefile"

        ( cd "$TRAL_EXT_SOFTWARE/$latestVer/core/" && make clean && make && make install ) # Installation
        rm -rf "$TRAL_EXT_SOFTWARE/"$latestVer""
    }
fi


# mafft is now accessible with the command "mafft"

######################
### Uninstalltion of MAFFT

# rm -rf "$TRAL_EXT_SOFTWARE"/mafft*

#      ( cd $INSTALLATION_PATH; \
# rm -f mafft linsi; rm -f mafft ginsi; rm -f mafft fftns; \
# rm -f mafft fftnsi; rm -f mafft nwns; rm -f mafft nwnsi; \
# rm -f mafft einsi; \
# rm -f mafft mafft-linsi; rm -f mafft mafft-ginsi; rm -f mafft mafft-fftns; \
# rm -f mafft mafft-fftnsi; rm -f mafft mafft-nwns; rm -f mafft mafft-nwnsi; \
# rm -f mafft mafft-einsi; rm -f mafft mafft-xinsi; rm -f mafft mafft-qinsi; \
# rm -f mafft mafft-distance; rm -f mafft mafft-profile; \
# rm -f mafft-homologs.rb; rm -f mafft-sparsecore.rb; \
# rm -rf libexec )
