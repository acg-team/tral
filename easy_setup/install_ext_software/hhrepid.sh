#!/bin/bash

# INSTALLING HHREPID #####
# HHrepID 1.0.0 – de novo Protein Repeat Detection
# HHrepID is a novel automated procedure for the de novo identification of repeats in protein sequences.
# It is able to detect the sequence signature of structural repeats in many proteins that have not yet been known to possess internal sequence symmetry, such as TIM barrels and outer membrane beta-barrels.
# HHrepID uses HMM-HMM comparison to exploit evolutionary information in the form of the multiple sequence alignment of homologs, but in contrast to HHrep, the new method has several novel characteristics:
# (1) automated generation of a multiple alignment of repeats;
# (2) utilization of the transitive nature of homology through a novel merging procedure based on the fully probabilistic treatment of alignments;
# (3) superior alignment quality through an algorithm that maximizes the expected accuracy of an alignment;
# (4) the ability to identify different repeats within complicated architectures or multiple domains through automatic domain boundary detection,
# (5) new statistical treatment yielding improved sensitivity.
# https://toolkit.tuebingen.mpg.de/#/tools/hhrepid

# Biegert A, Söding J:
# HHrepID: de novo protein repeat identification by probabilistic consistency
# Bioinformatics 2008 24(6), 807-814.

# Licence Summary:
# You are free
#   * to copy, distribute, display and perform the work
#   * to make derivative works

# Under the following conditions:
#   *  Attribution: You must give the original author credit.
#   *  Noncommercial: You may not use this work for commercial purposes.

# For any reuse or distribution, you must make clear to others the license terms of this work.

# Any of these conditions, in particular the second condition restricting the use for commercial purposes, 
# can be waived if you get permission from the copyright holder.

# Your fair use and other rights are in no way affected by the above.

# Full license can be found here: ftp://ftp.tuebingen.mpg.de/pub/protevo/HHrepID/LICENSE



######################
### Housekeeping

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P )
# other files are located one directory above
. $PARENT_PATH/configTRAL_path.cfg # provide paths from config file


######################
### Installation HHREPID

mkdir -p $TRAL_EXT_SOFTWARE/HHrepID



if [ ! -f $TRAL_EXT_SOFTWARE/HHrepID/hhrepid_64 ]; then # test if not already in directory
    LINK_HHREPID=ftp://ftp.tuebingen.mpg.de/pub/protevo/HHrepID # TRF Version 4.09 (Feb 22, 2016) Linux64
    wget $LINK_HHREPID/hhrepid_64 -P $TRAL_EXT_SOFTWARE/HHrepID   # download execution file
    wget $LINK_HHREPID/README -P $TRAL_EXT_SOFTWARE/HHrepID    # download README
fi

chmod +x $TRAL_EXT_SOFTWARE/HHrepID/hhrepid_64
cp $TRAL_EXT_SOFTWARE/HHrepID/hhrepid_64 /usr/local/bin/ # copy into system path
chmod +x /usr/local/bin/hhrepid_64 && echo -e  "\nHHrepID is in your system path /usr/local/bin/ and can be executed with the command hhrepid_64"

# HHrepID is now executable with hhrepid_64

######################
### Uninstall HHrepID

# rm -rf $TRAL_EXT_SOFTWARE/HHrepID/
# rm -rf /usr/local/bin/hhrepid_64