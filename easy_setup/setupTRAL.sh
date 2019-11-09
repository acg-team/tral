#!/usr/bin/env bash

# By execution of this script a little filesystem will be created within the INSTALLATION_PATH (default: /usr/local/bin).
# If you wish to change this path, do this within configTRAL_path.cfg
# Run this script as superuser if your INSTALLATION_PATH only can be accessed by root.


######################
# PREPARING FILESYSTEM AND INSTALLING TRAL
######################

shopt -s nocasematch # making comparisons case-insensitive
set -euo pipefail # exit on errors and undefined vars

if [[ ! "$1" =~ ^(setup|pip)$ ]]; then
    echo -e "\nPlease provide as argument how you want to install TRAL with setup.py (\"setup\" or \"pip\").\n"
    exit 1
fi

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg


# create needed directories to install tral
mkdir -p "$TRAL_PATH"
mkdir -p "$TRAL_EXT_SOFTWARE" # create directory for installation of external software

# directories will be added temporarely to PATH
[[ ":$PATH:" != *"$TRAL:$PATH"* ]] && PATH="$TRAL:$PATH"

# check for pip
if ! hash "${PIP3:-pip}"; then
    echo "${PIP3:-pip} is required. Set the PIP variable in configTRAL_path.cfg" >&2
    exit 1
fi


######################
### Installing TRAL


if [[ $1 == "setup" ]]; then

    ######################
    ### install tral with setup.py

    PARENT_PATH=$(  cd "$(dirname "${BASH_SOURCE[0]}")" || exit $?  ; cd .. ; pwd -P )

    # test if TRAL repository is already downloaded
    if [ -d "$PARENT_PATH/tral" ] && [ -e "$PARENT_PATH/setup.py" ] ; then
        echo -e "\nUsing TRAL source from $PARENT_PATH\n"

        (cd "$PARENT_PATH" && "${PYTHON3:-python3}" setup.py install ${TRAL_SETUP_ARGS:-}) || {
            echo -e "\nA problem occured while trying to install TRAL with \"${PYTHON3:-python3} setup.py install ${TRAL_SETUP_ARGS:-}\"."
            exit 1
        }
    else

        ### installing with git and setup.py
        # if the tral directory is not already downloaded, it has to be cloned from github
        echo -e "\nPlease clone the TRAL repository from github.\n"

        if [[ "${ACCEPT_ALL:-no}" = "yes" ]]; then
        yn=y
        else read -p "Do you want to clone the TRAL repository from github in "$TRAL_PATH"? yes(y) or no (n):" yn
        fi

        case $yn in
            [Yy]* )
                # check if git is installed
                if ! hash git 2>/dev/null; then
                    echo "Git is required to clone the repository of TRAL. Otherwise you can download the directory manually." >&2
                    exit 1
                fi

                git clone https://github.com/acg-team/tral.git "$TRAL_PATH/tral_repository" || {
                    echo -e "\nA problem occured while trying to clone the TRAL repository."
                    exit 1
                }
                (cd "$TRAL_PATH/tral_repository" && "${PYTHON3:-python3}" setup.py install ${TRAL_SETUP_ARGS:-})  || {
                    echo -e "\nA problem occured while trying to install TRAL with \"${PYTHON3:-python3} setup.py install ${TRAL_SETUP_ARGS:-}\"."
                    exit 1
                }
            ;;
            [Nn]* )
                echo -e "\nAbort."
                exit 1
            ;;
        esac
    fi



elif [[ $1 == "pip" ]]; then

    #####################
    ## install tral with pip
    #####################

    # -- TODO .tral will not automatically be added to $HOME

    {
        echo -e "start installation"
        "${PIP:-pip}" install tral|| {

        echo -e "\nA problem occured while trying to install TRAL with \"${PIP:-pip} install\"."
        exit $?
        }

    } && {

        echo -e "\n---------------------------------"
        echo -e "Installing TRAL with pip successful"
        echo -e "-----------------------------------\n"

    }


fi


######################
### Download p-Value distribution files

if [ ! -d "$TRAL_CONF/data/pvalue" ]; then
    echo -e  "\nIn order to calculate the p-Value of tandem repeat scores, available p-Value distributions need to be downloaded (2.6Gb) and placed in ~/.tral/data/pvalue."

    if [[ "${ACCEPT_ALL:-no}" = "yes" ]]; then
    yn=y
    else read -p "Would you like to do this? yes(y) or no (n):" yn
    fi

    case $yn in
        [Yy]* )
            {
                if [[ ! -f "$TRAL_CONF/data/pvalue.tar.gz" ]]; then
                    echo "DOWNLOADING"
                    wget "https://acg-team.ulozezoz.myhostpoint.ch/pvalue.tar.gz" -P "$TRAL_CONF/data"
                else
                    echo "SKIPPING DOWNLOAD"
                fi
                tar -xzf "$TRAL_CONF/data/pvalue.tar.gz" -C "$TRAL_CONF/data"
                rm -rf "$TRAL_CONF/data/pvalue.tar.gz"
                } || {
                exit 1
            }
        ;;
        [Nn]* )
            echo -e "\nYou can download this file later from ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper.\n"
        ;;
    esac
else
    echo -e "\nDirectory for p-Value distribution file already exists.\nWill not be updated."
fi

echo -e "\n---------------------------------"
echo -e "Installation of TRAL is completed."
echo -e "-----------------------------------\n"

######################
### Check if $TRAL_CONF was created during installation
if [ ! -d "$TRAL_CONF" ]; then
    echo -e ".tral directory was not created."
fi






