#!/bin/bash

# INSTALLING HMMER #####
# hmmer: Sequence profile model generation
# http://hmmer.org/documentation.html

# installing the newest version of hmmer with direct link
# HMMER will be obtained and compiled from source
# (!) HMMER will be installed within the $TRAL_EXT_SOFTWARE directory

######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file

[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE:$PATH"* ]] && PATH="$TRAL_EXT_SOFTWARE:$PATH"


######################
### Download and Installation hmmer

# grab the latest Version of hmmer for linux-intel-x86_64 (from the old ones)
sudo wget http://eddylab.org/software/hmmer/hmmer.tar.gz -P $TRAL_EXT_SOFTWARE 
sudo tar -xvzf $TRAL_EXT_SOFTWARE/hmmer.tar.gz -C $TRAL_EXT_SOFTWARE
latestVersion=$(tar tfz $TRAL_EXT_SOFTWARE/hmmer.tar.gz --exclude '*/*') # get name of version
cd $TRAL_EXT_SOFTWARE/$latestVersion

sudo ./configure --prefix $TRAL_EXT_SOFTWARE    # with this version, HMMER will be installed within the $TRAL_EXT_SOFTWARE directory
sudo make
sudo make check
sudo make install
cd ~
sudo rm -rf $TRAL_EXT_SOFTWARE/hmmer.tar.gz



