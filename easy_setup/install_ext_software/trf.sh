#!/usr/bin/env bash

# INSTALLING TRF #####
# TRF: Tandem repeat finder
# Tandem Repeats Finder is a program to locate and display tandem repeats in DNA sequences.
# In order to use the program, the user submits a sequence in FASTA format.
# https://tandem.bu.edu/trf/trf.unix.help.html

# G. Benson; "Tandem repeats finder: a program to analyze DNA sequences"; Nucleic Acids Research (1999); Vol. 27, No. 2, pp. 573-580.

# Licence:
# The author of this software grants to any individual or organization the right to use and to make an unlimited number of copies of this software.
# You may not de-compile, disassemble, reverse engineer, or modify the software.
# This software cannot be sold, incorporated into commercial software or redistributed.
# The author of this software accepts no responsibility for damages resulting from the use of this software and makes no warranty or representation, either express or implied, including but not limited to, any implied warranty of merchantability or fitness for a particular purpose.
# This software is provided as is, and the user assumes all risks when using it. 

######################
### Housekeeping

shopt -s nocasematch # making comparisons case-insensitive

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}


######################
### Installation TRF

if [ ! -f "$TRAL_EXT_SOFTWARE/trf409.linux64" ]; then # test if not already in directory
    LINK_TRF="tandem.bu.edu/trf/downloads/trf409.linux64" # TRF Version 4.09 (Feb 22, 2016) Linux64
    wget "$LINK_TRF" -P "$TRAL_EXT_SOFTWARE" || {   # download
    echo "Was not able to download TRF."
    exit $?
    }
else
    echo "TRF is already installed"
fi

chmod +x "$TRAL_EXT_SOFTWARE/trf409.linux64"
cp "$TRAL_EXT_SOFTWARE/trf409.linux64" "$INSTALLATION_PATH" # copy into system path
chmod +x "$INSTALLATION_PATH/trf409.linux64" 
{
if [ ! -h "$INSTALLATION_PATH/trf" ]; then    
    ln -s trf409.linux64 "$INSTALLATION_PATH/trf" # create symlink to executable file
fi
} && echo -e  "\nTRF is in your system path $INSTALLATION_PATH and can be executed with the command trf"

# TRF is now executable with trf409.linux64

######################
### Uninstall TRF

# rm -rf "$TRAL_EXT_SOFTWARE/trf409.linux64"
# rm -rf "$INSTALLATION_PATH/trf409.linux64"
# rm -rf "$INSTALLATION_PATH/trf"