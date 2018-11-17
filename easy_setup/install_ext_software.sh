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

    read -p "Would you like to install "$(basename "${software%%.sh}")"? Type \"y\" if YES:" y
    case $y in
        [Yy]* )
            echo ". install_ext_software/"$(basename "$software")""
            ;;
        * ) 
            echo -e "\nYou can install it later with the script $software.sh.\n"
            ;;
    esac

}

read -p "Would you like to install any external software? yes(y) or no (n):" yn
case $yn in
    [Yy]* )
        echo -e "\n"
        for software in install_ext_software/*.sh ; do install_ext_software $software ; done
        ;;
    [Nn]* ) 
        echo -e "\nNo external software will be installed right now."
        ;;
esac

