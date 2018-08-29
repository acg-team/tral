### T-REKS
source $HOME_USER/tral/bin/activate 
LINK_TREKS=http://bioinfo.montp.cnrs.fr/t-reks/T-Reks.jar
sudo wget $LINK_TREKS -P $TRAL_TOOLS
# Create an executable text file TREKS
echo '#!/bin/sh' >> TREKS
echo 'java -jar /usr/local/bin/tral_tools/T-Reks.jar "$@"' >> TREKS
chmod +x TREKS
sudo mv TREKS $TRAL_TOOLS/
# Add the TREKS file to the PATH variable
echo "export PATH=$PATH:$TRAL_TOOLS/TREKS" >> ~/.profile
source ~/.profile

