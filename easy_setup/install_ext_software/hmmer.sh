#!/usr/bin/env bash

# INSTALLING HMMER #####
# HMMER: Sequence profile model generation using profile hidden Markov models

# HMMER is used for searching sequence databases for sequence homologs, and for making sequence alignments.
# It implements methods using probabilistic models called profile hidden Markov models (profile HMMs).

# http://hmmer.org/documentation.html
# HMMER and its documentation are freely distributed under the #3-Clause BSD open source license.
# For a copy of the license, see opensource.org/licenses/BSD-3-Clause.


######################
### Housekeeping

shopt -s nocasematch # making comparisons case-insensitive
set -euo pipefail # exit on errors and undefined vars

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}

######################
### Download HMMER

## TODO: get full name of unzipped hmmer (with version)

for directory in "$TRAL_EXT_SOFTWARE/hmmer"*; do  # test if hmmer files not already in directory
    if [ -d "$directory" ]; then
        echo "HMMER is already downloaded and in "$directory""
    else
        {   
            LINK_HMMER="http://eddylab.org/software/hmmer/hmmer.tar.gz" # HHrepID Nov 22 2007
            wget "$LINK_HMMER" -P "$TRAL_EXT_SOFTWARE"   # download execution file
            tar -xvzf "$TRAL_EXT_SOFTWARE/hmmer.tar.gz" -C "$TRAL_EXT_SOFTWARE"
            rm -rf "$TRAL_EXT_SOFTWARE/hmmer.tar.gz"
            } || {
            echo "Couldn't download HMMER."
            exit $?
        }
        
    fi
done

######################
### Compile and Install HMMER

{
    (
        {
            cd "$TRAL_EXT_SOFTWARE/hmmer-"*
        } && {
            ./configure --prefix "$INSTALLATION_PATH"
            make clean
            make
            # "$INSTALLATION_PATH"/bin make check        # run a test suite
            make install
        } && {
            echo "Installation of HMMER done."
            ln -s "$INSTALLATION_PATH/bin/hmmbuild" "$INSTALLATION_PATH/hmmbuild"
            echo -e  "\nhmmbuild is in your path $INSTALLATION_PATH\n"
        }
    )
    } || {
    echo "Couldn't compile and install HMMER. Try to manually install the development version of HMMER by running hmmer-devel.sh which provides little endian support."
    exit $?
}

###############i#######
### Uninstall HMMER

## TODO: find other way to delete hmmer and adapt it within uninstall script!!

# {
#     ( cd "$TRAL_EXT_SOFTWARE/hmmer-"* && make uninstall )
#     rm -rf "$TRAL_EXT_SOFTWARE/hmmer-"* 
#     rm -rf "$INSTALLATION_PATH/hmmbuild"
# } || {
#     echo -e "HMMER removed."
# }
