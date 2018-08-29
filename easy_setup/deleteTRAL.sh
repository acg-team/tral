#!/bin/bash

#provide paths from config file (has to be in the same directory than setupTRAL.sh)
. configTRAL_path.cfg 

# TODO-- delete path variable
# TODO-- change all variable names
# TODO-- delete external software
# TODO-- 
# TODO-- 
# TODO-- 





# activate the python3 virtualenv
source $HOME/tral/bin/activate
# store all installed packages in requirements format in a text file
pip freeze > requirements.txt
# uninstall all the installed packages in the virtualenv
pip uninstall -r requirements.txt -y
# deactivate the virtualenv
deactivate
# remove the requirements.txt file
rm requirements.txt
echo'----------------------------------------------'
echo 'The python3 based virtualenv of tral is reset.'
echo'----------------------------------------------'

# same procedure as above but for the python2 virtualenv
source $HOME/tral2/bin/activate
pip freeze > requirements.txt
pip uninstall -r requirements.txt -y
deactivate
rm requirements.txt
echo'----------------------------------------------'
echo 'The python2 based virtualenv of tral is reset.'
echo'----------------------------------------------'

# delete the virtualenvironments
sudo rm -rf $HOME/tral
echo'----------------------------------------------'
echo 'The python3 based virtualenv is deleted.'
sudo rm -rf $HOME/tral2
echo'----------------------------------------------'
echo 'The python2 based virtualenv is deleted.'
echo'----------------------------------------------'

# delete all the config files of the virtualenvs
#sudo rm -rf $HOME/.tral
#sudo rm -rf $HOME/.tral2

# remove the variables set to the system PATH variable
directory_TRED=/usr/bin/tral_tools/TRED
directory_TREKS=/usr/bin/tral_tools/TREKS
directory_TRUST=/usr/bin/tral_tools/TRUST
PATH=:$PATH:
PATH=${PATH//:$directory_TRED:/:}
PATH=${PATH#:}; PATH=${PATH%:}
PATH=:$PATH:
PATH=${PATH//:$directory_TREKS:/:}
PATH=${PATH#:}; PATH=${PATH%:}
PATH=:$PATH:
PATH=${PATH//:$directory_TRUST:/:}
PATH=${PATH#:}; PATH=${PATH%:}
echo'----------------------------------------------'
echo 'Your system path variable is reset to:' $PATH
echo'----------------------------------------------'

# delete the directory with external software used by TRAL
sudo rm -rf /usr/bin/tral_tools
echo'----------------------------------------------'
echo 'All external software of TRAL is deleted.'
echo'----------------------------------------------'




echo 'Uninstallation of TRAL successfully completed.'
