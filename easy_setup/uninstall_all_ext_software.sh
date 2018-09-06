#!/bin/bash

# TODO -- add the remaining external software
# TODO -- add some comments

######################
### Explanation

# This script uninstalls all external software of TRAL for which a installation script is provided within easy_setup
# A proper uninstallation is only provided if the software is installed in the default path
# If you like to keep an installed program, outcomment the associated part in this script and run
# It does not uninstall TRAL or change the filesystem


######################
### Housekeeping

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 


######################
### Uninstallation of all external software

### ALF

rm -rf /usr/local/bin/darwin*
rm -rf /usr/local/bin/alfsim
rm -rf /usr/local/share/alfdarwin
rm -rf /usr/local/bin/alfdarwin*

### HMMER

apt-get remove hmmer

### MAFFT

dpkg --remove mafft
rm -rf $TRAL_EXT_SOFTWARE/mafft*

### PHOBOS

rm -rf $TRAL_EXT_SOFTWARE/phobos*
rm -rf /usr/local/bin/phobos*

### TREDPARSE

rm -rf $TRAL_ENV/python2
rm $TRAL_EXT_SOFTWARE/tred
rm /usr/local/bin/tred

### TREKS

rm -rf $TRAL_EXT_SOFTWARE/T-Reks.jar
rm -rf $TRAL_EXT_SOFTWARE/TREKS
rm -rf /usr/local/bin/TREKS

### TRF

rm -rf $TRAL_EXT_SOFTWARE/trf409.linux64
rm -rf /usr/local/bin/trf409.linux64

### TRUST

rm -rf $TRAL_EXT_SOFTWARE/TRUST_Align
rm -rf /usr/local/bin/TRUST

