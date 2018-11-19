#!/bin/bash


######################
# UNINSTALLING TRAL AND ITS EXTERNAL SOFTWARE
######################

######################
### Prepare Filesystem

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg


######################
### Delete TRAL

# delete tral path system
if [[ "$ACCEPT_ALL" = "yes" ]] || [[ "$ACCEPT_ALL" = "Yes" ]]; then
    yn=y
else read -p "Are you sure to delete TRAL? yes(y) or no (n):" yn
fi

case $yn in
    [Yy]* )
        rm -rf "$TRAL_PATH" || {
            echo -e "Was not able to delete TRAL"
            exit 1
        }
    ;;
    [Nn]* )
        echo -e ".tral was not removed."
    ;;
esac


# delete all the config files of the virtualenvs
if [[ "$ACCEPT_ALL" = "yes" ]] || [[ "$ACCEPT_ALL" = "Yes" ]]; then
    yn=y
else read -p "Do you wish to delete all configuration files and data as well? yes(y) or no (n):" yn
fi

case $yn in
    [Yy]* )
        rm -rf "$TRAL_CONF" || {
            echo -e "Was not able to delete .tral"
            exit 1
        }
    ;;
    [Nn]* )
        echo -e ".tral was not removed."
    ;;
esac

echo -e "---------------------------------------------"
echo -e "TRAL and its whole pathsystem is now deleted."
echo -e "---------------------------------------------\n"


######################
### Delete external software of TRAL

# delete the directory with external software used by TRAL
if [[ "$ACCEPT_ALL" = "yes" ]] || [[ "$ACCEPT_ALL" = "Yes" ]]; then
    yn=y
else read -p "Do you wish to uninstall all external software as well? yes(y) or no (n):" yn
fi

case $yn in
    [Yy]* )
        {
            . uninstall_all_ext_software.sh &&
            rm -rf "$TRAL_EXT_SOFTWARE"
            echo -e "\n----------------------------------------------"
            echo -e "External Software deleted."
            echo -e "----------------------------------------------\n"
            echo -e "\n----------------------------------------------"
            echo -e 'Uninstallation of TRAL successfully completed.'
            echo -e "----------------------------------------------\n"
            } || {
            echo -e "\nA problem occured while trying to delete the external software."
            exit 1
        }
    ;;
    [Nn]* )
        echo -e "\nExternal Software was not deleted and can (if it was available before) still be found in $TRAL_EXT_SOFTWARE."
        echo -e "\n----------------------------------------------"
        echo -e 'Uninstallation of TRAL successfully completed.'
        echo -e "----------------------------------------------\n"
    ;;
esac