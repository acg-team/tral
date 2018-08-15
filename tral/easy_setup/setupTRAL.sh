#!/bin/bash
######################
### Housekeeping
######################
sudo mkdir -p /usr/bin/tral_tools # -p create it only if it doesn't already exist
TRAL_TOOLS=/usr/bin/tral_tools

######################
### install virutalenv and activate
######################
# check if virtual env is installed
if [! hash virtualenv 2>/dev/null]; then

	echo >&2 "Install virtualenv by: sudo apt install virutalenv"
	exit 1
fi

# install virtual environment called "tral" with python3.5
virtualenv $HOME/tral -p python3.5
TRAL_PATH=$HOME/tral/
TRAL_CONF_PATH=$HOME/.tral/
mkdir -p $TRAL_PATH/output # -p create it only if it doesn't already exist 

# activate the virtual environment
source $TRAL_PATH/bin/activate 

######################
### install tral
######################
#$TRAL_PATH/bin/pip install tral
pip install tral

######################
### install ipython3 in virtualenv
######################

pip install ipython

######################
### install external software
######################
### hmmer: Sequence profile model generation
# download and installation
LINK_HMMER=http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
sudo wget $LINK_HMMER -P $TRAL_TOOLS
sudo tar zxf $TRAL_TOOLS/hmmer-3.1b2-linux-intel-x86_64.tar.gz -C $TRAL_TOOLS 
$TRAL_TOOLS/hmmer-3.1b2-linux-intel-x86_64/configure
make
make check
sudo rm -rf $TRAL_TOOLS/hmmer-3.1b2-linux-intel-x86_64.tar.gz
# export PATH=foo:$PATH @TODO: change foo here! Check where file is. Wrapper necessary?

### mafft: Alignment of tandem repeat units
LINK_MAFFT=https://mafft.cbrc.jp/alignment/software/mafft_7.397-1_amd64.deb
sudo wget $LINK_MAFFT -P $TRAL_TOOLS
sudo dpkg -i $TRAL_TOOLS/mafft_7.397-1_amd64.deb
sudo rm -rf $TRAL_TOOLS/mafft_7.397-1_amd64.deb

### PHOBOS
LINK_PHOBOS=http://www.rub.de/ecoevo/cm/phobos-v3.3.12-linux.tar.gz
sudo wget $LINK_PHOBOS -P $TRAL_TOOLS
sudo tar zxf $TRAL_TOOLS/phobos-v3.3.12-linux.tar.gz -C $TRAL_TOOLS 
sudo rm -rf $TRAL_TOOLS/phobos-v3.3.12-linux.tar.gz
#run by:
#/usr/bin/tral_tools/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 
# export PATH=foo:$PATH @TODO: change foo here! Check where file is. Wrapper necessary?

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

### T-REKS
source $HOME_USER/tral/bin/activate 
LINK_TREKS=http://bioinfo.montp.cnrs.fr/t-reks/T-Reks.jar
sudo wget $LINK_TREKS -P $TRAL_TOOLS
# Create an executable text file TREKS
echo '#!/bin/sh' >> TREKS
echo 'java -jar /usr/bin/tral_tools/T-Reks.jar "$@"' >> TREKS
chmod +x TREKS
sudo mv TREKS $TRAL_TOOLS/
# Add the TREKS file to the PATH variable
echo "export PATH=$PATH:$TRAL_TOOLS/TREKS" >> ~/.profile
source ~/.profile

### TRF
#download manually from: https://tandem.bu.edu/trf/trf409.linux64.download.html
sudo mv $HOME/Downloads/trf409.linux64 $TRAL_TOOLS/
chmod +x trf409.linux64

### TRUST
LINK_TRUST=http://www.ibi.vu.nl/programs/trustwww/trust.tgz
sudo wget $LINK_TRUST -P $TRAL_TOOLS
sudo tar -xvzf $TRAL_TOOLS/trust.tgz -C $TRAL_TOOLS 
# Create an executable text file TRUST
echo '#!/bin/sh' >> TRUST
echo 'java -Xmx30G -cp /usr/bin/tral_tools/Align nl.vu.cs.align.SelfSimilarity "$@"' >> TRUST
chmod +x TRUST
sudo mv TRUST $TRAL_TOOLS/
# Add the TRUST file to the PATH variable
echo "export PATH=$PATH:$TRAL_TOOLS/TRUST" >> ~/.profile
source ~/.profile

### ALF
LINK_ALF=http://abacus.gene.ucl.ac.uk/daniel/alf/ALF_standalone.tar.gz
sudo wget $LINK_ALF -P $TRAL_TOOLS
sudo tar -xvzf $TRAL_TOOLS/ALF_standalone.tar.gz -C $TRAL_TOOLS 
sudo rm -rf $TRAL_TOOLS/ALF_standalone.tar.gz
# By default, this will install the binary and script files in /usr/local/bin and the lib directory in /usr/local/share/alfdarwin:
cd $TRAL_TOOLS/ALF_standalone/
sudo ./install.sh
cd ~
# ALF is executable with alfdarwin

######################
### Download p-Value distribution files
######################
sudo wget ftp://ftp.vital-it.ch/papers/vital-it/Bioinformatics-Schaper/pvalue.tar.gz -P $TRAL_CONF_PATH/data/
sudo tar -xzf pvalue.tar.gz -C $TRAL_CONF_PATH/data/
sudo rm -rf $TRAL_CONF_PATH/data/pvalue.tar.gz
