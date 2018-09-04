#!/bin/bash

# here comes a nice description how to use setupTRAL.sh and other shell scripts in this directory

# run this script as a superuser (depending on where TRAL and depencies should be installed)

######################
# PREPARING FILESYSTEM AND INSTALLING TRAL
######################

# TODO-- how to use within script (bash requires sudo)
# TODO-- adapt the other setupfiles to the same changes!
# TODO-- install TRAL with pip (tral should be installed anyway)
# TODO-- describe how to use setupTRAL.sh and scripts for external software, and deleteTRAL.sh
# TODO-- add installation path for external software and install external software within this folder
# TODO-- create README for easy_setup (or include in main README of TRAL?)


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
### install virtualenv and activate

# check if virtualenv is installed
while [ hash virtualenv 2>/dev/null ] ; do
    echo "Installing virtualenv required by: apt-get install virtualenv."
    read -p "Do you wish to install this program? yes(y) or no (n):" yn
        case $yn in
            [Yy]* )
                apt-get install virtualenv || { 
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

# create virtual environment called "python3_5" with python3.5
virtualenv $TRAL_ENV/python3_5 -p python3.5

# activate the virtual environment
. $TRAL_ENV/python3_5/bin/activate

######################
### install ipython3 in virtualenv
pip install ipython

######################
### installing TRAL

# TODO -- to use this script, the TRAL repository should already be downloaded (?)

if [[ $1 == "git" ]]; then

    ######################
    ### install tral with git and setup.py
    
    # check if git is installed
    while [ hash git 2>/dev/null ] ; do
        echo "Git is required to install TRAL."
    done
    
    # installing with git
    echo -e "\nPlease clone the TRAL repository from github into $TRAL_PATH.\n"
    read -p "Do you want to clone the TRAL repository from github in the given directory? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            git clone https://github.com/acg-team/tral.git $TRAL
            (cd $TRAL && python $TRAL/setup.py develop)                   ## TODO change from develop to install
            ;;
        [Nn]* ) 
            echo -e "\nAbort."
            exit 1
            ;;
    esac

    
elif [[ $1 == "pip" ]]; then

    #####################
    ## install tral with pip
    #####################

    # -- TODO .tral will not automatically be added to $HOME

    echo -e "start installation"
    $TRAL_ENV/python3_5/bin/pip3 install --target=$TRAL_ENV/python3_5/bin/ tral && 

    echo -e "\n---------------------------------"
    echo -e "Installing TRAL with pip successful"
    echo -e "-----------------------------------\n"

else
    echo -e "\nPlease give as argument if you want to install TRAL with git or with pip."
    exit 1


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
            echo -e "\nYou can download this files later from ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper"
            ;;
    esac
else
    echo -e "\nDirectory for p-Value distribution files already exists.\nWill not be updated."
fi


