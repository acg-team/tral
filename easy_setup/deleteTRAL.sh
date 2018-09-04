#!/bin/bash

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 

# TODO-- delete unneeded parts
# TODO-- comment
# TODO-- include error handling
# TODO-- make style consistent

# TODO-- can we not just delete the folder $TRAL_ENV? therefore i deleted all the part with freezing virtenv


# delete tral path system
rm -rf $TRAL_PATH

# delete all the config files of the virtualenvs
rm -rf $TRAL_CONF
#rm -rf $HOME/.tral2

echo -e "---------------------------------------------"
echo -e "TRAL and its whole pathsystem is now deleted."
echo -e "---------------------------------------------\n"


# delete the directory with external software used by TRAL
read -p "Do you wish to uninstall all external software as well? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            rm -rf $TRAL_EXT_SOFTWARE && echo -e "\nExternal Software deleted." || { 
            echo -e "\nA problem occured while trying to delete the external software."
            exit 1 
            }
            ;;
        [Nn]* ) 
            echo -e "\nExternal Software was not deleted and can (if it was available before) still be found in $TRAL_EXT_SOFTWARE."
            ;;
    esac

echo -e "\n----------------------------------------------"
echo -e 'Uninstallation of TRAL successfully completed.'
echo -e "----------------------------------------------\n"