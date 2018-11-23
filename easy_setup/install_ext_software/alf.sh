#!/bin/bash

# INSTALLING ALF #####

# ALF—A Simulation Framework for Genome Evolution

# ALF simulates a root genome into a number of related genomes.
# Result files include the resulting gene sequences, true tree and true MSAs.
# A description of ALF can be found in the following article:
# Daniel A Dalquen, Maria Anisimova, Gaston H Gonnet, Christophe Dessimoz: ALF - A Simulation Framework for Genome Evolution. Mol Biol Evol, 29(4):1115-1123, April 2012. http://mbe.oxfordjournals.org/content/29/4/1115
# If you use ALF for your publication, please cite.

######################
### Housekeeping

shopt -s nocasematch # making comparisons case-insensitive

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}

[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE/bin:$PATH"* ]] && export PATH="$TRAL_EXT_SOFTWARE/bin:$PATH"


######################
### Download and Installation ALF

if [ ! -d "$TRAL_EXT_SOFTWARE/ALF_standalone" ]; then # test if not already in directory
    {
    LINK_ALF="http://abacus.gene.ucl.ac.uk/daniel/alf/ALF_standalone.tar.gz"
    wget "$LINK_ALF" -P "$TRAL_EXT_SOFTWARE"    # download
    tar -xvzf "$TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz" -C "$TRAL_EXT_SOFTWARE"
    } || {
        echo "Couldn't download or unzip ALF."
        exit $?
    }

fi

rm -rf "$TRAL_EXT_SOFTWARE/ALF_standalone.tar.gz"
(cd "$TRAL_EXT_SOFTWARE/ALF_standalone" && "$TRAL_EXT_SOFTWARE/ALF_standalone/install.sh") # installation of ALF


######################
### Uninstall ALF (default paths!)

# rm -rf /usr/local/bin/alfdarwin.linux64
# rm -rf /usr/local/bin/alfsim
# rm -rf /usr/local/share/alfdarwin


