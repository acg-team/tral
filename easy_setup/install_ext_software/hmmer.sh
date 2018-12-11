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

shopt -s nocasematch # making comparisons case-insensitive

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}


APT_CMD=$(which apt)
YUM_CMD=$(which yum)
DNF_CMD=$(which dnf)
CONDA_CMD=$(which conda)
BREW_CMD=$(which brew)
PORT_CMD=$(which port)

######################
### Installation HMMER
# installing the newest version of hmmer depending on your system


if [[ ! -z $APT_CMD ]]; then            # Linux (Ubuntu, Debian...)
    apt install hmmer
elif [[ ! -z $DNF_CMD ]]; then          # Linux (Fedora)
    dnf install hmmer
elif [[ ! -z $YUM_CMD ]]; then          # Linux (older Fedora)
    yum install hmmer
elif [[ ! -z $CONDA_CMD ]]; then        # Anaconda
    conda install -c bioconda hmmer
elif [[ ! -z $BREW_CMD ]]; then         # OS/X, HomeBrew
    brew install hmmer
elif [[ ! -z $PORT_CMD ]]; then         # OS/X, MacPorts
    port install hmmer
else
    echo "Error: Unfortunatelly it is not possible to install HMMER on your system."
    exit 1;
fi

######################
### Uninstall HMMER

# if [[ ! -z $APT_CMD ]]; then            # Linux (Ubuntu, Debian...)
#     apt remove hmmer
# elif [[ ! -z $DNF_CMD ]]; then          # Linux (Fedora)
#     dnf remove hmmer
# elif [[ ! -z $YUM_CMD ]]; then          # Linux (older Fedora)
#     yum remove hmmer
# elif [[ ! -z $CONDA_CMD ]]; then        # Anaconda
#     conda remove -c bioconda hmmer
# elif [[ ! -z $BREW_CMD ]]; then         # OS/X, HomeBrew
#     brew remove hmmer
# elif [[ ! -z $PORT_CMD ]]; then         # OS/X, MacPorts
#     port uninstall hmmer
# else
#     echo "Error: Something went wrong while trying to remove HMMER, maybe it is already uninstalled."
#     exit 1;
# fi