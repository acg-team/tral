#!/bin/bash

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
rm -rf "$TRAL_EXT_SOFTWARE"/mafft*

### PHOBOS

rm -rf "$TRAL_EXT_SOFTWARE"/phobos*
rm -rf "$INSTALLATION_PATH"/phobos*

### TREDPARSE

rm -rf "$TRAL_ENV/python2"
rm "$TRAL_EXT_SOFTWARE/tred"
rm "$INSTALLATION_PATH/tred"

### TREKS

rm -rf "$TRAL_EXT_SOFTWARE/T-Reks.jar"
rm -rf "$TRAL_EXT_SOFTWARE/T-REKS"
rm -rf "$INSTALLATION_PATH/T-REKS"

### TRF

rm -rf "$TRAL_EXT_SOFTWARE/trf409.linux64"
rm -rf "$INSTALLATION_PATH/trf409.linux64"

### TRUST

rm -rf "$TRAL_EXT_SOFTWARE/TRUST_Align"
rm -rf "$INSTALLATION_PATH/TRUST"

### HHrepID

rm -rf "$TRAL_EXT_SOFTWARE/HHrepID/"
rm -rf "$INSTALLATION_PATH/hhrepid_64"


### XSTREAM

rm -rf "$TRAL_EXT_SOFTWARE/XSTREAM"
rm -rf "$INSTALLATION_PATH/XSTREAM"

