# easy_setup for TRAL

Attention:
This README is just a first draft as well as the whole installation procedure. Adaption is coming soon.

These scripts will help you to easily install TRAL and its dependencies without going to deep into detail about the different installation procedures.
It automatically sets up a little filesystem, where a virtual environment (virtenv) for TRAL will be ininitalized.
Moreover, you can decide which external software should be downloaded and installed without the need to read through their (sometimes complicated) installation procedures.

## Getting Started

To get started with TRAL, simply download the directory easy_setup and execute the bash scripts as described below.
You can adapt the default installing path within configTRAL_path.cfg. Please do only some changes to the variable INSTALLATION_PATH.

### Prerequisites

Until easy_setup only works for linux64. If you use another operating system, you can have a closer look into the scripts and adapt them for your specific OS.
To install TRAL and its external software within usr/local/bin (default) you need to execute the scripts as root.
In case you want to install TRAL by automatically download the repository you need to have installed git.


### Installation of TRAL

To install TRAL within its virtual environment simply run the script setupTRAL.sh with either the argument "pip" or "git".

## Installation of TRAL with "git"

At time the version provide in "pip" is not todate, therefore it is recommended to install TRAL with "git".

```
sudo ./setupTRAL.sh
```

Unless you have git installed, you may need to download the github repository manually (https://github.com/acg-team/tral.git).


### Installation of external software for TRAL

After setting up TRAL you can install external software. The script setupTRAL.sh automatically ask you wheter you want to install any of the recommended software.
By pressing y/n you can decide wheter you want to install it now or later.

If you decide to install the external software later, you can run a script for each individual external software which you can find within the directory setup_tral/install_ext_software.
It automatically will be installed within the directory INSTALLATION_PATH/tral_external_software.

At the moment scripts for installing these external software are available:

- alf
- hmmer
- mafft
- phobos
- tredparse
- t-reks
- trf
- trust

In the end of each installation script for the external software you can find how to uninstall the software.
Either outcomment all installation part and activate the uninstallation part or run the commands directly in the commandline.


### Activation and Use of TRAL

Since TRAL is installed within a virtual environment, you first have to activate it.
For that you can run the script activateTRAL either with

```
. activateTRAL.sh
```
or
```
source activateTRAL.sh
```.


### Uninstallation of TRAL and its external software

To uninstall TRAL (and is external software if you wish) run the script deleteTRAL.sh.



## Author

Paulina Naef



