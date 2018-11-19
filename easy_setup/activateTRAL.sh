#!/bin/bash

# script to run whenever one want to start TRAL
# run this script with ". activateTRAL.sh" or "source activateTRAL.sh"
# to deactivate the virtual tral environment type "deactivate" into the shell


# provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg

# activate the virtual environment
. "$TRAL_ENV/python3/bin/activate"  || {
    echo -e "\nA problem occured while trying to activate the virtual environment."
    exit 1
}

# add directory with tral and with external software to path (beginning)

[[ ":$PATH:" != *"$TRAL:$PATH"* ]] && export PATH="$TRAL:$PATH"
[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE:$PATH"* ]] && export PATH="$TRAL_EXT_SOFTWARE:$PATH"
[[ ":$PATH:" != *"$TRAL_EXT_SOFTWARE/bin:$PATH"* ]] && export PATH="$TRAL_EXT_SOFTWARE/bin:$PATH"
