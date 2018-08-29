#!/bin/bash

# INSTALLING HMMER #####
# hmmer: Sequence profile model generation
# http://hmmer.org/documentation.html


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Download and Installation hmmer

# grab the latest Version of hmmer for linux-intel-x86_64 (from the old ones)
latestVer=$(wget -qO- http://eddylab.org/software/hmmer3/CURRENT/ |
            grep linux-intel-x86_64 | 
            sed -n 's/.*href="\([^"]*\).*/\1/p')
sudo wget http://eddylab.org/software/hmmer3/CURRENT/${latestVer} -P $TRAL_EXT
sudo tar -xvzf $TRAL_EXT/${latestVer} -C $TRAL_EXT 
cd $TRAL_EXT/${latestVer%.tar.gz}
sudo ./configure --prefix $TRAL_SOFTWARE
sudo make
sudo make check
sudo make install
cd ~
sudo rm -rf $TRAL_EXT/${latestVer}
#export PATH=foo:$PATH @TODO: change foo here! Check where file is. Wrapper necessary?
