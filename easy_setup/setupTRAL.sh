#!/bin/bash

# By execution of this script a little filesystem will be created within the INSTALLATION_PATH (default: /usr/local/bin).
# If you wish to change this path, do this within configTRAL_path.cfg
# TRAL will be installed within a virtual environment (virtenv).
# Run this script as as superuser if your INSTALLATION_PATH only can be accessed by root.


######################
# PREPARING FILESYSTEM AND INSTALLING TRAL
######################

if [[ ! "$1" =~ ^(setup|pip)$ ]]; then
    echo -e "\nPlease provide as argument how you want to install TRAL with setup.py \"setup\" or \"pip\").\n"
    exit 1
fi

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 

# create needed directories to install tral
mkdir -p $TRAL_PATH/{tral,tral_env,output}
mkdir -p $TRAL_EXT_SOFTWARE # create directory for installation of external software

# directories will be added temporarely to PATH 
[[ ":$PATH:" != *"$TRAL:$PATH"* ]] && PATH="$TRAL:$PATH"

######################
### Install virtualenv and activate

# check for pip
if ! hash "${PIP3:-pip}"; then
    echo "${PIP3:-pip} is required. Set the PIP variable in configTRAL_path.cfg" >&2
    exit 1
fi

# check if virtualenv is installed
while [ ! -x $(which virtualenv 2>/dev/null) ]; do
    echo "Installing virtualenv required by: "${PIP:-pip}" install virtualenv."
    read -p "Do you wish to install this program? yes(y) or no (n):" yn
        case $yn in
            [Yy]* )
                "${PIP3:-pip}" install virtualenv || { 
                echo -e "\nA problem occured while trying to install virtualenv."
                exit 1 
                }
                ;;
            [Nn]* ) 
                echo -e "\nAbort."
                exit 1
                ;;
        esac
done

# create virtual environment called "python3" with python3.5
virtualenv $TRAL_ENV/python3 -p $PYTHON3 || exit $?

# activate the virtual environment
. $TRAL_ENV/python3/bin/activate || exit $?

######################
### Installing TRAL




if [[ $1 == "setup" ]]; then

    ######################
    ### install tral with setup.py
    
    PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
    CURRENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P ) # do not use only pwd (shortcuts could have other path)

    # test if TRAL repository is already downloaded
    if [[ "${CURRENT_PATH%%/easy_setup}" == $PARENT_PATH ]] && [ -e "$PARENT_PATH/setup.spy" ] ; then
        python $PARENT_PATH/setup.py install
    else

        ### installing with git and setup.py
        # if the tral directory is not already downloaded, it has to be cloned from github
        echo -e "\nPlease clone the TRAL repository from github.\n"
        read -p "Do you want to clone the TRAL repository from github in the home directory? yes(y) or no (n):" yn
        case $yn in
            [Yy]* )
                # check if git is installed
                if ! hash git 2>/dev/null; then
                    echo "Git is required to clone the repository of TRAL. Otherwise you can download the directory manually." >&2
                    exit 1
                fi

                git clone https://github.com/acg-team/tral.git $HOME
                (cd $HOME && python tral/setup.py install)
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

    echo -e "start installation"
    $TRAL_ENV/python3/bin/pip3 install --target=$TRAL_ENV/python3/bin/ tral && 

    echo -e "\n---------------------------------"
    echo -e "Installing TRAL with pip successful"
    echo -e "-----------------------------------\n"


fi

deactivate


######################
### Download p-Value distribution files

if [ ! -d "$TRAL_CONF/data/pvalue" ]; then
    echo -e  "\nIn order to calculate the p-Value of tandem repeat scores, available p-Value distributions need to be downloaded (2.6Gb) and placed in ~/.tral/data/pvalue."
    read -p "Would you like to do this? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            wget ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper/pvalue.tar.gz -P $TRAL_CONF/data
            tar -xzf $TRAL_CONF/data/pvalue.tar.gz -C $TRAL_CONF/data
            rm -rf $TRAL_CONF/data/pvalue.tar.gz
            ;;
        [Nn]* ) 
            echo -e "\nYou can download this files later from ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper.\n"
            ;;
    esac
else
    echo -e "\nDirectory for p-Value distribution files already exists.\nWill not be updated."
fi

######################
### Installation of external software for TRAL


. install_ext_software.sh # script to install external software within the same directory will be executed

echo -e "\n---------------------------------"
echo -e "Installation of TRAL is completed."
echo -e "-----------------------------------\n"






