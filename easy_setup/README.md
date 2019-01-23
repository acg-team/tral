# easy_setup for TRAL

These scripts will help you to easily install TRAL and its dependencies without going to deep into detail of the different installation procedures.
It automatically sets up a little filesystem, where a virtual environment (virtenv) for TRAL will be initalized.
Moreover, you can decide which external software should be downloaded and installed without the need to read through their (sometimes complicated) installation procedures.

## Getting Started

To get started with TRAL, simply clone the github directory https://github.com/acg-team/tral.git or download the directory easy_setup and then execute the bash scripts as described below.
You can adapt the default installing path within configTRAL_path.cfg. Please only change the variables FILES and INSTALLATION_PATH.

## Prerequisites

These installation scripts can be used on all Linux systems, but still require python3, pip, unzip and git.
Per default TRAL and its external software will be installed in usr/local/bin. Then you may execute the scripts as root.
If you wish to change this path you can adjust $INSTALLATION_PATH in configTRAL.cfg

These setup scripts require python 3, pip, unzip, and git and run on any UNIX based operation system.

## Installation of TRAL

To install TRAL within you can simply run the script setupTRAL.sh with either the argument "pip" or "git".
It is recommended to install TRAL within a virtual environment (e. g. virtualenv).

### Installation of TRAL with "git"

Currently the version provided in "pypi" is not up-to-date, therefore it is recommended to install TRAL with python setup.py install "setup".

```
sudo ./setupTRAL.sh setup
```
If you haven't git installed, you may need install git or to download the github repository manually (https://github.com/acg-team/tral.git).


### Installation of external software for TRAL

After installing TRAL with setupTRAL.sh you can install external software. 
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
Either comment out all installation part and uncomment uninstallation part or run the commands directly in the commandline.
 
You may need to run these scripts as root ($INSTALLATION_PATH is usr/local/bin per default).


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
