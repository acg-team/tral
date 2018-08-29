### PHOBOS
LINK_PHOBOS=http://www.rub.de/ecoevo/cm/phobos-v3.3.12-linux.tar.gz
sudo wget $LINK_PHOBOS -P $TRAL_TOOLS
sudo tar zxf $TRAL_TOOLS/phobos-v3.3.12-linux.tar.gz -C $TRAL_TOOLS 
sudo rm -rf $TRAL_TOOLS/phobos-v3.3.12-linux.tar.gz
#run by:
#/usr/local/bin/tral_tools/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 
# export PATH=foo:$PATH @TODO: change foo here! Check where file is. Wrapper necessary?
