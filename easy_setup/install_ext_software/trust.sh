### TRUST
LINK_TRUST=http://www.ibi.vu.nl/programs/trustwww/trust.tgz
sudo wget $LINK_TRUST -P $TRAL_TOOLS
sudo tar -xvzf $TRAL_TOOLS/trust.tgz -C $TRAL_TOOLS 
# Create an executable text file TRUST
echo '#!/bin/sh' >> TRUST
echo 'java -Xmx30G -cp /usr/local/bin/tral_tools/Align nl.vu.cs.align.SelfSimilarity "$@"' >> TRUST
chmod +x TRUST
sudo mv TRUST $TRAL_TOOLS/
# Add the TRUST file to the PATH variable
echo "export PATH=$PATH:$TRAL_TOOLS/TRUST" >> ~/.profile
source ~/.profile
