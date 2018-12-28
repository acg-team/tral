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
shopt -s nocasematch # making comparisons case-insensitive


######################
### Uninstallation of all external software

### ALF

rm -rf /usr/local/bin//alfdarwin.linux64
rm -rf /usr/local/bin/alfsim
rm -rf /usr/local/share/alfdarwin


### MAFFT

rm -rf "$TRAL_EXT_SOFTWARE"/mafft*

     ( cd $INSTALLATION_PATH; \
rm -f mafft linsi; rm -f mafft ginsi; rm -f mafft fftns; \
rm -f mafft fftnsi; rm -f mafft nwns; rm -f mafft nwnsi; \
rm -f mafft einsi; \
rm -f mafft mafft-linsi; rm -f mafft mafft-ginsi; rm -f mafft mafft-fftns; \
rm -f mafft mafft-fftnsi; rm -f mafft mafft-nwns; rm -f mafft mafft-nwnsi; \
rm -f mafft mafft-einsi; rm -f mafft mafft-xinsi; rm -f mafft mafft-qinsi; \
rm -f mafft mafft-distance; rm -f mafft mafft-profile; \
rm -f mafft-homologs.rb; rm -f mafft-sparsecore.rb; \
rm -rf libexec )

### PHOBOS

rm -rf "$TRAL_EXT_SOFTWARE/phobos_64_libstdc++6"
rm -rf "$INSTALLATION_PATH/phobos_64_libstdc++6"
rm -rf "$INSTALLATION_PATH/phobos"

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
rm -rf "$INSTALLATION_PATH/trf"

### TRUST

rm -rf "$TRAL_EXT_SOFTWARE/TRUST_Align"
rm -rf "$INSTALLATION_PATH/TRUST"

### HHrepID

rm -rf "$TRAL_EXT_SOFTWARE/HHrepID/"
rm -rf "$INSTALLATION_PATH/hhrepid_64"

### XSTREAM

rm -rf "$TRAL_EXT_SOFTWARE/XSTREAM"
rm -rf "$INSTALLATION_PATH/XSTREAM"

### HMMER

if [[ ! -z $APT_CMD ]]; then            # Linux (Ubuntu, Debian...)
    apt remove hmmer
elif [[ ! -z $DNF_CMD ]]; then          # Linux (Fedora)
    dnf remove hmmer
elif [[ ! -z $YUM_CMD ]]; then          # Linux (older Fedora)
    yum remove hmmer
elif [[ ! -z $CONDA_CMD ]]; then        # Anaconda
    conda remove bioconda hmmer
elif [[ ! -z $BREW_CMD ]]; then         # OS/X, HomeBrew
    brew remove hmmer
elif [[ ! -z $PORT_CMD ]]; then         # OS/X, MacPorts
    port uninstall hmmer
else
    echo "Error: Something went wrong while trying to remove HMMER, maybe it is already uninstalled."
    exit 1;
fi

