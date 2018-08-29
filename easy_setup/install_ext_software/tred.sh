### TRED
#only supported for python 2.7, python 3 not yet supported
# install virtual environment called "tral2" with python2.7
virtualenv $HOME/tral2 -p python2.7
# activate the virtual environment
source $HOME/tral2/bin/activate
sudo apt install ipython
# install tral
$HOME/tral2/bin/pip install tral
# install TRED
pip install tredparse
# Wrapper, sourcing tral2, call tred.py with all arguments "$@"
# Create an executable text file TRED
echo '#!/bin/sh' >> TRED
echo 'source $HOME/tral2/bin/activate' >> TRED
echo 'tred.py "$@"' >> TRED
chmod +x TRED
sudo mv TRED $TRAL_TOOLS/
# Add the TRED file to the PATH variable
echo "export PATH=$PATH:$TRAL_TOOLS/TRED" >> ~/.profile
source ~/.profile
deactivate
