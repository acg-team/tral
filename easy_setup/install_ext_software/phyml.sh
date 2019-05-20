#!/bin/bash

# INSTALLING PhyML #####
# PhyML is a software package which primary task that is to estimate maximum
# likelihood phylogenies from alignments of nucleotide or amino acid sequences.

# Copyright 1999 - 2008 by PhyML Development Team.
# The software PhyML is provided “as is” without warranty of any kind. In no event shall
# the authors or his employer be held responsible for any damage resulting from the use
# of this software, including but not limited to the frustration that you may experience in
# using the package. All parts of the source and documentation except where indicated are
# distributed under the GNU public licence. See http://www.opensource.org for details


# PhyML: “A simple, fast and accurate algorithm to estimate large phylogenies
# by maximum likelihood” Guindon S., Gascuel O. 2003, Systematic Biology,
# 52(5):696-704


# "New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0."
# Guindon S., Dufayard J.F., Lefort V., Anisimova M., Hordijk W., Gascuel O.
# Systematic Biology, 59(3):307-21, 2010.


######################
### Housekeeping
 
PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" # provide paths from config file


######################
### Installation PhyML


if [ ! -d "$TRAL_EXT_SOFTWARE/PhyML-3.1" ]; then # test if not already in directory
    LINK_PhyML="http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip"
    wget "$LINK_PhyML" -P "$TRAL_EXT_SOFTWARE"
    unzip "$TRAL_EXT_SOFTWARE/PhyML-3.1.zip" -d "$TRAL_EXT_SOFTWARE"
    rm -rf "$TRAL_EXT_SOFTWARE/PhyML-3.1.zip"
fi


# Create an executable file PhyML

cp "$TRAL_EXT_SOFTWARE/PhyML-3.1/PhyML-3.1_linux64" "$INSTALLATION_PATH"  # copy wrapper file to execute PhyML into system path
chmod +x "$INSTALLATION_PATH/PhyML-3.1_linux64" && echo -e "\nPhyML is in your system path "$INSTALLATION_PATH" and can be executed with the command \"PhyML-3.1_linux64\""

# PhyML is executable with the command PhyML-3.1_linux64


######################
### Uninstall PhyML

# rm -rf "$TRAL_EXT_SOFTWARE/PhyML-3.1"
# rm -rf "$INSTALLATION_PATH/PhyML-3.1_linux64"