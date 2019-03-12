#!/usr/bin/env bash

# The user should have reading/writing rights on the system folders (i.e. /usr/local, /usr/local/include).

# INSTALLING Castor #####

######################
### Housekeeping

shopt -s nocasematch # making comparisons case-insensitive

PARENT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; cd .. ; pwd -P ) # other files are located one directory above
. "$PARENT_PATH/configTRAL_path.cfg" || {  # provide paths from config file
    echo "configTRAL_path.cfg not found"
    exit $?
}

mkdir -p "$TRAL_EXT_SOFTWARE/castor"

######################
### Compiling and installing the dependencies and castor

(
    cd "$TRAL_EXT_SOFTWARE/castor"
    {       
        {
            ## bpp-core http://biopp.univ-montp2.fr/
            (
            git clone https://github.com/BioPP/bpp-core
            cd bpp-core
            git checkout tags/v2.4.0 -b v240
            mkdir build
            cd build
            cmake ..
            # if one does not want to install with default path:
            # cmake -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" ..
            make clean
            sudo make install
            )
        } || { 
            echo "A problem occured while trying to install bpp-core."
            exit $? 
        }
    } && {
        {
            ## bpp-seq http://biopp.univ-montp2.fr/
            (
            git clone https://github.com/BioPP/bpp-seq
            cd bpp-seq
            git checkout tags/v2.4.0 -b v240
            mkdir build
            cd build
            cmake ..
            # if one does not want to install with default path:
            # cmake -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" ..
            make clean
            sudo make install
            )
        } || { 
            echo "A problem occured while trying to install bpp-seq."
            exit $? 
        }
    } && {
        {
            ## bpp-phy http://biopp.univ-montp2.fr/
            (
            git clone https://github.com/BioPP/bpp-phyl
            cd bpp-phyl
            git checkout tags/v2.4.0 -b v240
            mkdir build
            cd build
            cmake ..
            # if one does not want to install with default path:
            # cmake -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" ..
            make clean
            sudo make install
            )
        } || { 
            echo "A problem occured while trying to install bpp-phy."
            exit $? 
        }
    } && {
        {
            ## boost - C++ Libraries http://www.boost.org/
            (
            wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
            tar xvf boost_1_66_0.tar.gz
            rm -r boost_1_66_0.tar.gz
            cd boost_1_66_0
            ./bootstrap.sh
            ./b2
            sudo ./b2 install
            )
        } || { 
            echo "A problem occured while trying to install boost."
            exit $? 
        }
    } && {
        {  
            ## glog - Google Logging Library https://github.com/google/glog
            (
            git clone https://github.com/google/glog
            cd glog
            cmake -H. -Bbuild -G "Unix Makefiles"
            sudo cmake --build build --target install
            )
        } || { 
            echo "A problem occured while trying to install glog."
            exit $? 
        }
    } && {  
        {
            ## TSHLib - Tree Search Heuristics Library 
            (
            git clone https://github.com/acg-team/tshlib.git # do not use the bitbucket
            # username for github required!
            cd tshlib
            git checkout develop
            cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
            sudo make install
            )
        } || { 
            echo "A problem occured while trying to install TSHLib."
            exit $? 
        }
    } && {
        {
            ## Intel TBB - Intel(R) Threading Building Blocks 2018
            (
            cd /opt
            sudo wget https://github.com/01org/tbb/releases/download/2018_U5/tbb2018_20180618oss_lin.tgz
            sudo tar -xvf tbb2018_20180618oss_lin.tgz
            )
        } || { 
            echo "A problem occured while trying to install Intel TBB."
            exit $? 
        }
    } && {
        {
            ## compiling Castor with dynamic linking
            (
            git clone https://github.com/acg-team/castor.git
            cd castor
            cmake --target castor -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
            sudo make
            )
        } || { 
            echo "A problem occured while trying to compile Castor."
            exit $? 
        }
    }
)

