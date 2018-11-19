#!/bin/bash

# INSTALLING TRUST #####
# TRUST: Tracking repeats using significance and transitivity.
# A method for ab-initio determination of internal repeats in proteins.
# The high sensitivity and accuracy of the method is achieved by exploiting the concept of transitivity of alignments.
# ISMB/ECCB 2004 conference (Glasgow, UK), appeared in Bioinformatics. 2004 Aug 4;20 Suppl 1:i311-i317. 
# http://www.ibi.vu.nl/programs/trustwww/


######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}

######################
### Installation TRUST

mkdir -p "$TRAL_EXT_SOFTWARE/TRUST_Align"
if [ ! -d "$TRAL_EXT_SOFTWARE/TRUST_Align/Align" ]; then # test if not already in directory
    LINK_TRUST="http://www.ibi.vu.nl/programs/trustwww/trust.tgz"
    wget "$LINK_TRUST" -P "$TRAL_EXT_SOFTWARE/TRUST_Align"
    tar -xvzf "$TRAL_EXT_SOFTWARE/TRUST_Align/trust.tgz" -C "$TRAL_EXT_SOFTWARE/TRUST_Align"
    rm -rf "$TRAL_EXT_SOFTWARE/TRUST_Align/trust.tgz" || {
        echo "Was not able to download or unzip TRUST"
        exit $?
    }
fi



# Create an executable file TRUST
Align="$TRAL_EXT_SOFTWARE/TRUST_Align/Align"

echo '#!/bin/sh
# wrapper file to easily start TRUST

java -Xmx30G -cp' $Align ' nl.vu.cs.align.SelfSimilarity "$@"' > "$TRAL_EXT_SOFTWARE/TRUST_Align/TRUST"
chmod +x "$TRAL_EXT_SOFTWARE/TRUST_Align/TRUST"
cp "$TRAL_EXT_SOFTWARE/TRUST_Align/TRUST" "$INSTALLATION_PATH"  # copy wrapper file to execute TRUST into system path
chmod +x "$INSTALLATION_PATH/TRUST" && echo -e "\nTRUST is in your system path $INSTALLATION_PATH and can be executed with the command \"TRUST\""

# TRUST is executable with the command TRUST


######################
### Uninstall TRUST

# rm -rf "$TRAL_EXT_SOFTWARE/TRUST_Align"
# rm -rf "$INSTALLATION_PATH/TRUST"