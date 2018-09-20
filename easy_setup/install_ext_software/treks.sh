#!/bin/bash

# INSTALLING T-REKS #####
# T-REKS: identification of Tandem Repeats in sequences with a K-meanS based algorithm. 
# T-REKS is an algorithm for de novo detection and alignment of repeats in sequences based on K-means algorithm.
# Minimal length of repeat arrays is 9 for true homorepeats and 14 for other repeats with potential biological meaning.
# T-reks: identification of tandem repeats in sequences with a k-means based algorithm. Jorda J, Kajava AV(2009). Bioinformatics 25 (20), 2632-2638 (https://www.ncbi.nlm.nih.gov/pubmed/19671691)


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file



######################
### Installation T-REKS

LINK_TREKS=https://dali.crbm.cnrs.fr/tools/treks/T-Reks.jar
wget --no-check-certificate $LINK_TREKS -P $TRAL_EXT_SOFTWARE # no check certificate because otherwise the download will fail

if [ ! -f $TRAL_EXT_SOFTWARE/T-Reks.jar ]; then # test if T-Reks.jar is downloaded or manually put into the directory.
    echo "Something went wrong with downloading T-REKS"

else # Create an executable file for TREKS
    echo '#!/bin/sh
    # wrapper file to easily start T-REKS

    java -jar' $TRAL_EXT_SOFTWARE/T-Reks.jar ' "$@"' > $TRAL_EXT_SOFTWARE/TREKS
    chmod +x $TRAL_EXT_SOFTWARE/TREKS
    cp $TRAL_EXT_SOFTWARE/TREKS /usr/local/bin/  # copy wrapper file to execute T-REKS into system path
    chmod +x /usr/local/bin/TREKS && echo -e "\nTREKS is in your system path /usr/local/bin/ and can be executed with the command \"TREKS\""


# TREKS is executalbe with the command "TREKS"
fi

######################
### Uninstallation of T-REKS

# rm -rf $TRAL_EXT_SOFTWARE/T-Reks.jar
# rm -rf $TRAL_EXT_SOFTWARE/TREKS
# rm -rf /usr/local/bin/TREKS

