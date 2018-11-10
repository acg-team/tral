#!/bin/bash

# Run this script after setting up TRAL by using setupTRAL.sh.
# This script will automatically called within setupTRAL or can be executed afterwards.
# You may need superuser privileges.

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 


######################
### Installing external software


install_ext_software () {
    for var in "$@"
    do
        read -p "Would you like to install $var? Type \"y\" if YES:" y
        case $y in
            [Yy]* )
                . install_ext_software/$var.sh
                ;;
            * ) 
                echo -e "\nYou can install it later with the script $var.sh.\n"
                ;;
        esac
    done
}

read -p "Would you like to install any external software? yes(y) or no (n):" yn
case $yn in
    [Yy]* )
        echo -e "\n"
        install_ext_software alf hhrepid hmmer mafft phobos tredparse treks trf trust xstream phyml     
        ;;
    [Nn]* ) 
        echo -e "\nNo external software will be installed right now."
        ;;
esac
