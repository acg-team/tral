#!/bin/bash

# INSTALLING HMMER #####
# HMMER: Sequence profile model generation using profile hidden Markov models

# HMMER is used for searching sequence databases for sequence homologs, and for making sequence alignments.
# It implements methods using probabilistic models called profile hidden Markov models (profile HMMs).

# http://hmmer.org/documentation.html
# HMMER and its documentation are freely distributed under the #3-Clause BSD open source license.
# For a copy of the license, see opensource.org/licenses/BSD-3-Clause.


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Installation HMMER

# installing the newest version of hmmer with a direct link

apt-get update
apt-get upgrade
apt install hmmer

# to install  a specific version number
# apt install hmmer=<version_number> 


######################
### Uninstall HMMER

# apt-get remove hmmer