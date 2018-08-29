#!/bin/bash

######################
# PREPARING FILESYSTEM AND INSTALLING TRAL
######################

# TODO-- how to use sudo within script (bash requires sudo)
# TODO-- adapt the other setupfiles to the same changes!
# TODO-- install TRAL with pip
# TODO-- describe how to use setupTRAL.sh and scripts for external software, and deleteTRAL.sh
# TODO-- add installation path for external software and install external software within this folder


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
    echo "Installing virtualenv required by: sudo apt-get install virtualenv."
    read -p "Do you wish to install this program? yes(y) or no (n):" yn
        case $yn in
            [Yy]* )
                sudo apt-get install virtualenv || { 
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

# create virtual environment called "tral_env" with python3.5
virtualenv $TRAL_ENV -p python3.5

# activate the virtual environment
source $TRAL_ENV/bin/activate

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

    
# else
#     #####################
#     ## install tral with pip
#     #####################
#     $TRAL_PATH/bin/pip install tral


#         check if pip is installed
#        while [ hash pip 2>/dev/null ] ; do
#            echo "Pip is required to install TRAL "
#        done

#        sudo pip install tral || { echo -e "\nA problem occured while trying to install TRAL." ; exit 1 }


fi

######################
### Download p-Value distribution files

if [ ! -d "$TRAL_CONF/data/pvalue" ]; then
    echo -e  "\nIn order to calculate the p-Value of tandem repeat scores, available p-Value distributions need to be downloaded (2.6Gb) and placed in ~/.tral/data/pvalue."
    read -p "Would you like to do this? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            sudo wget ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper/pvalue.tar.gz -P $TRAL_CONF/data
            sudo tar -xzf $TRAL_CONF/data/pvalue.tar.gz -C $TRAL_CONF/data
            sudo rm -rf $TRAL_CONF/data/pvalue.tar.gz
            ;;
        [Nn]* ) 
            echo -e "\nYou can download this files later from ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper"
            ;;
    esac
else
    echo -e "\nDirectory for p-Value distribution files already exists.\nWill not be updated."
fi


