#!/bin/bash

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 

# TODO-- delete unneeded parts
# TODO-- comment
# TODO-- include error handling
# TODO-- make style consistent


# activate the python3 virtualenv
source $TRAL_ENV/bin/activate
# store all installed packages in requirements format in a text file
sudo pip freeze > requirements.txt
# uninstall all the installed packages in the virtualenv
sudo pip uninstall -r requirements.txt -y
# deactivate the virtualenv
deactivate
# remove the requirements.txt file
sudo rm requirements.txt
echo '----------------------------------------------'
echo 'The python3 based virtualenv of tral is reset.'
echo '----------------------------------------------'

# # TODO -- decide if tral environment for python2 necessary
# source $HOME/tral2/bin/activate
# pip freeze > requirements.txt
# pip uninstall -r requirements.txt -y
# deactivate
# rm requirements.txt
# echo'----------------------------------------------'
# echo 'The python2 based virtualenv of tral is reset.'
# echo'----------------------------------------------'

# delete tral path system
sudo rm -rf $TRAL_PATH

# delete the virtualenvironments
# sudo rm -rf $TRAL_ENV
# echo '----------------------------------------------'
# echo 'The python3 based virtualenv is deleted.'
# sudo rm -rf $HOME/tral2
# echo '----------------------------------------------'
# echo 'The python2 based virtualenv is deleted.'
# echo'---------------------------------------------'

# delete all the config files of the virtualenvs
sudo rm -rf $TRAL_CONF
#sudo rm -rf $HOME/.tral2






#### better do not system wide changes to the path variable (does this work at all??)
# # remove the variables set to the system PATH variable
# directory_TRED=/usr/bin/tral_tools/TRED
# directory_TREKS=/usr/bin/tral_tools/TREKS
# directory_TRUST=/usr/bin/tral_tools/TRUST
# PATH=:$PATH:
# PATH=${PATH//:$directory_TRED:/:}
# PATH=${PATH#:}; PATH=${PATH%:}
# PATH=:$PATH:
# PATH=${PATH//:$directory_TREKS:/:}
# PATH=${PATH#:}; PATH=${PATH%:}
# PATH=:$PATH:
# PATH=${PATH//:$directory_TRUST:/:}
# PATH=${PATH#:}; PATH=${PATH%:}
# echo'----------------------------------------------'
# echo 'Your system path variable is reset to:' $PATH
# echo'----------------------------------------------'

# delete the directory with external software used by TRAL
read -p "Do you wish to uninstall all external software as well? yes(y) or no (n):" yn
    case $yn in
        [Yy]* )
            sudo rm -rf $TRAL_EXT_SOFTWARE && echo -e "\nExternal Software deleted." || { 
            echo -e "\nA problem occured while trying to delete the external software."
            exit 1 
            }
            ;;
        [Nn]* ) 
            echo -e "\nExternal Software was not deleted and can (if not already deleted) still be found in $TRAL_EXT_SOFTWARE."
            exit 1
            ;;
    esac

echo -e '\nUninstallation of TRAL successfully completed.'
