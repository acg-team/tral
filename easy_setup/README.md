# easy_setup for TRAL

These scripts will help you to easily install TRAL and its dependencies without going to deep into detail of the different installation procedures.
It automatically sets up a little filesystem, where a virtual environment (virtenv) for TRAL will be initalized.
Moreover, you can decide which external software should be downloaded and installed without the need to read through their (sometimes complicated) installation procedures.

## Getting Started

To get started with TRAL, simply clone the github directory https://github.com/acg-team/tral.git or download the directory easy_setup and then execute the bash scripts as described below.
You can adapt the default installing path within configTRAL_path.cfg. Please only change the variables FILES and INSTALLATION_PATH.

## Prerequisites

These setup scripts require python 3, pip, unzip, and git and run on any UNIX based operation system.

To install TRAL and its external software within usr/local/bin (default for variable $INSTALLATION_PATH) you need to execute the scripts as root.
In case you want to install TRAL from the git repository automatically, you need to have git installed on your computer.

## Installation of TRAL

To install TRAL within its virtual environment you can simply run the script setupTRAL.sh with either the argument "pip" or "git".

### Installation of TRAL with "git"

Currently the version provided in "pypi" is not up-to-date, therefore it is recommended to install TRAL with python setup.py install "setup".

```
sudo ./setupTRAL.sh setup
```
If you haven't git installed, you may need install git or to download the github repository manually (https://github.com/acg-team/tral.git).


### Installation of external software for TRAL

After setting up TRAL with setupTRAL.sh you can install external software. 
The script install_ext_software.sh automatically iterates through all installation scripts for each recommended external software of TRAL and ask you for each single software if you want to install it (press y/n). Please confirm to their respective licence.

Otherwise, you can run an installation script for each external software individually which can be found within the directory setup_tral/install_ext_software.
The software will be downloaded in the directory $FILES/tral_external_software.
Currently, installation scripts for the following external software are available:

- alf
- hmmer
- hhrepid
- mafft
- phobos
- tredparse  (virtualenv required to run in python2)
- t-reks
- trf
- trust
- xstream
- tred

On the bottom of each installation script for the external software an uninstallation procedure can be found.
Either comment all installation part and uncomment uninstallation part or run the commands directly in the commandline.

It may be that you need to be root to install external software.

## Activation and Use of TRAL

It is recommended to install TRAL within a virtual environment (e. g. virtualenv).

## Uninstallation of TRAL and its external software

To uninstall TRAL (and is external software if you wish) run the script deleteTRAL.sh.

```
\.deleteTRAL.sh
```
You will be asked if you want to uninstall the external software as well.

To uninstall only all installed external software run:

```
\.uninstall_all_ext_software.sh
```
