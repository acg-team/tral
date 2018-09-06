#!/bin/bash

# INSTALLING T-REKS #####
# T-REKS: identification of Tandem REpeats in sequences with a K-meanS based algorithm. (https://www.ncbi.nlm.nih.gov/pubmed/19671691)
# 
# TODO -- Website (http://bioinfo.montp.cnrs.fr/) is down, any possibility to download it somewhere else? 
#       Mailed to: julien.jorda@crbm.cnrs.fr & andrey.kajava@crbm.cnrs.fr
#       Unable to reach Juline Jorda
#       Tried to get in contact with CRBM (http://www.crbm.cnrs.fr/en/contact/)

# Dear Sir or Madam

# I'm a student at the ZHAW in Switzerland (Applied Computational Life Sciences) and would like to use your Tandem Repeat Detector T-REKS (published 2009) for my Master Thesis.

# Unfortunately, I see no more option to download your program since your (old?) website (http://bioinfo.montp.cnrs.fr) is down.
# The programe was available with http://bioinfo.montp.cnrs.fr/?r=t-reks.
# Is there any public version available somewhere else? I wasn't able to find one. Or is the program no longer maintained?

# Many thanks in advance.


# Kind regards,
# Paulina Näf



######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file



######################
### Installation T-REKS

# ATTTENTION: Not working since the website for downloading is down
# LINK_TREKS=http://bioinfo.montp.cnrs.fr/t-reks/T-Reks.jar
# sudo wget $LINK_TREKS -P $TRAL_EXT_SOFTWARE

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
### Uninstall T-REKS

# rm -rf $TRAL_EXT_SOFTWARE/T-Reks.jar
# rm -rf $TRAL_EXT_SOFTWARE/TREKS
# rm -rf /usr/local/bin/TREKS

