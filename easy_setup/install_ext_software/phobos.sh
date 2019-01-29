#!/usr/bin/env bash

# INSTALLING PHOBOS #####
# Phobos - a tandem repeat search tool for complete genomes
# http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm
# Academic (non-commercial) version of Phobos, Copyright (c) Christoph Mayer 2006-2017

# License:
# This program Phobos is copyright protected. It is only distributed from the authors web page (www.rub.de/spezzoo/cm/cm_phobos.htm).
# For academic and non-commercial usage this program can be used free of charge. Results obtained
# with this program can be published without restrictions, provided the program and its author are acknowledged by name.
# A commercial license for Phobos can be obtained from the author.

# TODO -- change default way to access PHOBOS in the config.ini file (phobos_*)

######################
### Housekeeping

shopt -s nocasematch # making comparisons case-insensitive

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}

######################
### Download and Installation PHOBOS

## TODO: Does the test work if phobos not installed in a system path?

if [ ! -x "$(command -v phobos_64_libstdc++6)" ]; then # test if not already in directory
    
    if [[ "$ACCEPT_ALL" = "yes" ]]; then
        ma=a
    else
        {
            echo -e  "\nThis program Phobos is copyright protected. It is only distributed from the authors web page (www.rub.de/spezzoo/cm/cm_phobos.htm).\n"
            echo -e  "\nWould you like to register and download the program manually? A popup will open."
            echo -e  "Otherwise you can directly download the phobos-v3.3.12-linux version automatically.\n"
            read -p "manually (m) or automatically (a)?:" ma
        }
    fi
    case $ma in
        [Mm]* )                                                                 # redirects user to download page
            if which xdg-open > /dev/null
            then
                xdg-open "https://www.ruhr-uni-bochum.de/ecoevo/cm/regist_form.htm"
            elif which gnome-open > /dev/null
            then
                gnome-open "https://www.ruhr-uni-bochum.de/ecoevo/cm/regist_form.htm"
            fi
            
            echo -e  "\nAfter downloading your version of choice you can unzip it and put the binaries into your PATH.\n"
        ;;
        [Aa]* )
            echo -e "\nDo you confirm to the authors copyright (Copyright (c) Christoph Mayer 2006-2017)?\n This program is for academic and non-commercial usage only."
            
            if [[ "$ACCEPT_ALL" = "yes" ]]; then
                yn=y
            else read -p "Would you like to download the phobos-v3.3.12-linux version? yes(y) or no (n):" yn
            fi
            
            case $yn in
                [Yy]* )
                    {
                    LINK_PHOBOS="http://www.rub.de/ecoevo/cm/phobos-v3.3.12-linux.tar.gz"
                    wget "$LINK_PHOBOS" -P "$TRAL_EXT_SOFTWARE"
                    tar zxf "$TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux.tar.gz" -C "$TRAL_EXT_SOFTWARE"
                    rm -rf "$TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux.tar.gz"
                    cp "$TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6" "$INSTALLATION_PATH" # copies binaries to system path $INSTALLATION_PATH/
                    if [ ! -h "$INSTALLATION_PATH/phobos" ]; then 
                        ln -s "$INSTALLATION_PATH/phobos_64_libstdc++6" "$INSTALLATION_PATH/phobos" # create symlink to executable file
                    fi
                    } || {
                        echo -e "Something went wrong with downloading or installing phobos."
                        exit $?
                    }
                    echo -e "\nPHOBOS binaries are now in the user path $INSTALLATION_PATH"
                    ;;
                [Nn]* )
                    echo -e "Abort."
                ;;
            esac
        ;;
    esac
    
else
    echo "PHOBOS already installed." && exit 0
    
fi



# Per default phobos is accessible via phobos (as long as the binary is copied to the system path).


######################
### Uninstall PHOBOS

# rm -rf "$TRAL_EXT_SOFTWARE/phobos_64_libstdc++6"
# rm -rf "$INSTALLATION_PATH/phobos_64_libstdc++6"
# rm -rf "$INSTALLATION_PATH/phobos"


