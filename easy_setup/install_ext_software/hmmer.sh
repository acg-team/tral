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
### Installation HMMER

sudo apt-get update
sudo apt-get upgrade
sudo apt install hmmer

# to install  a specific version number
# sudo apt install hmmer=<version_number> 


######################
### Uninstall HMMER

# sudo apt-get remove hmmer