#!/bin/bash

# INSTALLING HMMER #####
# hmmer: Sequence profile model generation
# http://hmmer.org/documentation.html

# installing the newest version of hmmer with direct link


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Download and Installation hmmer

# grab the latest Version of hmmer for linux-intel-x86_64 (from the old ones)
sudo wget http://eddylab.org/software/hmmer/hmmer.tar.gz -P $TRAL_EXT
sudo tar -xvzf $TRAL_EXT/hmmer.tar.gz -C $TRAL_EXT
latestVersion=$(tar tfz $TRAL_EXT/hmmer.tar.gz --exclude '*/*')
cd $TRAL_EXT/$latestVersion
sudo ./configure --prefix $TRAL_SOFTWARE
sudo make
sudo make check
sudo make install
cd ~
sudo rm -rf $TRAL_EXT/hmmer.tar.gz
#export PATH=foo:$PATH @TODO: change foo here! Check where file is. Wrapper necessary?

