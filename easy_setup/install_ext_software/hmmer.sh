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

apt-get update
apt-get upgrade
apt install hmmer

# to install  a specific version number
# apt install hmmer=<version_number> 


######################
### Uninstall HMMER

# apt-get remove hmmer