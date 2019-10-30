#!/usr/bin/env bash

# INSTALLING XSTREAM #####
# XSTREAM: variable sequence tandem repeats extraction and architecture modeling
# XSTREAM is a tool for rapidly identifying and modeling the architecture of “fundamental” Tandem Repeats (TRs) in protein sequences.
# Due to the general nature of TRs, however, any sequence including DNA (or even numbers!) can be processed. 
# XSTREAM is freely available only for academic non-profit use and may not be distributed or copied without permission.
# For questions regarding XSTREAM availability, please contact the XSTREAM team.

# Newman, Aaron & B Cooper, James. (2007). XSTREAM: A practical algorithm for identification and architecture modeling of tandem repeats in protein sequences. BMC bioinformatics. 8. 382. 10.1186/1471-2105-8-382. 
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-382


######################
### Housekeeping

# Check if Java is installed
if java -version 2>&1 >/dev/null | grep "java version\|openjdk version" ; then   
	echo "Java is already installed."; 
else   
	echo "Java NOT installed!"
	exit $?;
fi



LINK_XSTREAM="https://amnewmanlab.stanford.edu/xstream/XSTREAMprog/xstream.zip"

shopt -s nocasematch # making comparisons case-insensitive
set -euo pipefail # exit on errors and undefined vars

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}


######################
### Installation XSTREAM

mkdir -p "$TRAL_EXT_SOFTWARE/XSTREAM"
if [ ! -f "$TRAL_EXT_SOFTWARE/XSTREAM/xstream.jar" ]; then # test if not already in directory
    wget "$LINK_XSTREAM" -P "$TRAL_EXT_SOFTWARE/XSTREAM"
    unzip "$TRAL_EXT_SOFTWARE/XSTREAM/xstream.zip" -d "$TRAL_EXT_SOFTWARE/XSTREAM"
    rm -rf "$TRAL_EXT_SOFTWARE/XSTREAM/xstream.zip" || {
        echo "Was not able to download or unzip XSTREAM"
        exit $?
    }
fi


# Create an executable file XSTREAM

echo '#!/bin/sh
# wrapper file to easily start XSTREAM

java -jar' "$TRAL_EXT_SOFTWARE/XSTREAM/xstream.jar" ' "$@"' > "$TRAL_EXT_SOFTWARE/XSTREAM/XSTREAM"
chmod +x "$TRAL_EXT_SOFTWARE/XSTREAM/XSTREAM"
cp "$TRAL_EXT_SOFTWARE/XSTREAM/XSTREAM" "$INSTALLATION_PATH"  # copy wrapper file to execute XSTREAM into system path
chmod +x "$INSTALLATION_PATH/XSTREAM" && echo -e "\nXSTREAM is in your path $INSTALLATION_PATH and can be executed with the command \"XSTREAM\""

# XSTREAM is executable with the command XSTREAM


######################
### Uninstall XSTREAM

# rm -rf "$TRAL_EXT_SOFTWARE/XSTREAM"
# rm -rf "$INSTALLATION_PATH/XSTREAM"
