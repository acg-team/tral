# easy_setup for TRAL

> **ATTENTION: This easy_setup is a first version which can have some unknown bugs. Moreover, the README may not be complete. Adaption is coming soon.**

These scripts will help you to easily install TRAL and its dependencies without going to deep into detail of the different installation procedures.
It automatically sets up a little filesystem, where a virtual environment (virtenv) for TRAL will be ininitalized.
Moreover, you can decide which external software should be downloaded and installed without the need to read through their (sometimes complicated) installation procedures.

## Getting Started

To get started with TRAL, simply download the directory easy_setup and execute the bash scripts as described below.
You can adapt the default installing path within configTRAL_path.cfg. Please only change the variable INSTALLATION_PATH.

## Prerequisites

easy_setup only works for linux64 by now. If you use another operating system, you can have a closer look into the scripts and adapt them for your specific OS.
To install TRAL and its external software within usr/local/bin (default) you need to execute the scripts as root.
In case you want to install TRAL from the git repository automatically, you need to have git installed on your computer.


## Installation of TRAL

To install TRAL within its virtual environment you can simply run the script setupTRAL.sh with either the argument "pip" or "git".

### Installation of TRAL with "git"

Currently the version provided in "pypi" is not up-to-date, therefore it is recommended to install TRAL with "git".

```
sudo ./setupTRAL.sh git
```
If you haven't git installed, you may need to download the github repository manually (https://github.com/acg-team/tral.git).


### Installation of external software for TRAL

After setting up TRAL with setupTRAL.sh you can install external software. The script setupTRAL.sh automatically ask you if you want to install any of the recommended software.
By pressing y/n you can decide wheter you want to install it now or later.

If you decide to install the external software later,you can install them later by executing the script install_ext_software. Otherwise, you can run an installation script for each external software individually which can be found within the directory setup_tral/install_ext_software.
The software will be installed within the directory INSTALLATION_PATH/tral_external_software.
Currently, installation scripts for the following external software are available:

- alf
- hmmer
- hhrepid
- mafft
- phobos
- tredparse
- t-reks
- trf
- trust

In the end of each installation script for the external software an uninstallation procedure can be found.
Either outcomment all installation part and uncomment uninstallation part or run the commands directly in the commandline.


## Activation and Use of TRAL

Since TRAL will be installed within a virtual environment, you need to activate it, each time you want to use TRAL.
For tactivation you can run the script activateTRAL.sh either with

```
. activateTRAL.sh
```
or
```
source activateTRAL.sh
```

## Uninstallation of TRAL and its external software

To uninstall TRAL (and is external software if you wish) run the script deleteTRAL.sh.



## Author

Paulina Naef

