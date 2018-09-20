#!/bin/bash


######################
# UNINSTALLING TRAL AND ITS EXTERNAL SOFTWARE
######################

# TODO-- include error handling

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 


######################
### Delete TRAL

# delete tral path system
rm -rf $TRAL_PATH

# delete all the config files of the virtualenvs
rm -rf $TRAL_CONF
#rm -rf $HOME/.tral2

echo -e "---------------------------------------------"
echo -e "TRAL and its whole pathsystem is now deleted."
echo -e "---------------------------------------------\n"


######################
### Delete external software of TRAL

# delete the directory with external software used by TRAL
read -p "Do you wish to uninstall all external software as well? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            . uninstall_all_ext_software.sh &&
            rm -rf $TRAL_EXT_SOFTWARE
            echo -e "\n----------------------------------------------" 
            echo -e "External Software deleted." 
            echo -e "----------------------------------------------\n"
            echo -e "\n----------------------------------------------"
            echo -e 'Uninstallation of TRAL successfully completed.'
            echo -e "----------------------------------------------\n" ||
            echo -e "\nA problem occured while trying to delete the external software."
            exit 1
            ;;
        [Nn]* ) 
            echo -e "\nExternal Software was not deleted and can (if it was available before) still be found in $TRAL_EXT_SOFTWARE."
            echo -e "\n----------------------------------------------"
            echo -e 'Uninstallation of TRAL successfully completed.'
            echo -e "----------------------------------------------\n" 
            ;;
    esac



