#!/bin/bash

# INSTALLING PHOBOS #####
# Phobos - a tandem repeat search tool for complete genomes
# http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm
# Download academic (non-commercial) version of Phobos, Copyright (c) Christoph Mayer 2006-2017

# License:
# This program Phobos is copyright protected. It is only distributed from the authors web page (www.rub.de/spezzoo/cm/cm_phobos.htm).
# For academic and non-commercial usage this program can be used free of charge. Results obtained
# with this program can be published without restrictions, provided the program and its author are acknowledged by name.
# A commercial license for Phobos can be obtained from the author.

# TODO -- PHOBOS should only be downloaded from the authors webpage. Therefore a popup will show up, which tells the user to download the distribution
# into their $TRAL_EXT_SOFTWARE
# TODO -- this script is only to help install, the authors right have to be respected


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Download and Installion PHOBOS

echo -e  "\nThis program Phobos is copyright protected. It is only distributed from the authors web page (www.rub.de/spezzoo/cm/cm_phobos.htm).\n"
echo -e  "\nWould you like to register and download the program manually? A popup will open."
echo -e  "Otherwise you can directly download the phobos-v3.3.12-linux version automatically."
read -p "manually (m) or automatically (a)?:" ma
case $ma in
    [Mm]* )                                                                 # redirects user to download page
        if which xdg-open > /dev/null
        then
            xdg-open https://www.ruhr-uni-bochum.de/ecoevo/cm/regist_form.htm
        elif which gnome-open > /dev/null
        then
            gnome-open https://www.ruhr-uni-bochum.de/ecoevo/cm/regist_form.htm
        fi

        echo -e  "\nAfter downloading your version of choice you can unzip it and put the binaries into your PATH.\n"
        ;;
    [Aa]* ) 
        echo -e  "\n"
        read -p "Do you confirm to the authors copyright and download the phobos-v3.3.12-linux version? yes(y) or no (n):" yn
        case $yn in
            [Yy]* )                                                            
                LINK_PHOBOS=http://www.rub.de/ecoevo/cm/phobos-v3.3.12-linux.tar.gz 
                sudo wget $LINK_PHOBOS -P $TRAL_EXT_SOFTWARE
                sudo tar zxf $TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux.tar.gz -C $TRAL_EXT_SOFTWARE
                sudo rm -rf $TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux.tar.gz
                sudo cp -n $TRAL_EXT_SOFTWARE/phobos-v3.3.12-linux/bin/phobos* /usr/local/bin/      # copies binaries to user path /usr/local/bin/
                echo -e "\nPHOBOS binaries are now in the user path /usr/local/bin/"
                ;;
            [Nn]* ) 
                echo -e "Abort."
                ;;
        esac
        ;;
esac



######################
### Uninstall PHOBOS

# sudo rm -rf $TRAL_EXT_SOFTWARE/phobos*
# sudo rm -rf /usr/local/bin/phobos*