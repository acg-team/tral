#!/usr/bin/env bash

# Run this script after setting up TRAL by using setupTRAL.sh.
# This script will automatically called within setupTRAL or can be executed afterwards.
# You may need superuser privileges.

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg
shopt -s nocasematch # making comparisons case-insensitive

######################
### Installing external software


install_ext_software () {
    
    if [[ "$ACCEPT_ALL" = "yes" ]]; then
        y=y
    else read -p "Would you like to install "$(basename "${software%%.sh}")"? Type \"y\" if YES:" y
    fi
    
    case $y in
        [Yy]* )
            . install_ext_software/"$(basename "$software")" || {
                echo -e "\nA problem occured while trying to run install_ext_software/"$(basename "$software")"."
                exit 1
            }
        ;;
        * )
            echo -e "\nYou can install it later with the script $software.sh.\n"
        ;;
    esac
    
}

if [[ "$ACCEPT_ALL" = "yes" ]]; then
    yn=y
else read -p "Would you like to install any external software? yes(y) or no (n):" yn
fi

case $yn in
    [Yy]* )
        echo -e "\n"
        for software in "install_ext_software"/*.sh
        do install_ext_software $software  || {
            echo -e "\nA problem occured while trying to call install_ext_software."
            exit 1
        }
        done
    ;;
    [Nn]* )
        echo -e "\nNo external software will be installed right now."
    ;;
esac

