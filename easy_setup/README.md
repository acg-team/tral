# easy_setup for TRAL

These scripts will help you to easily install TRAL and its dependencies without going too deep into detail of the different installation procedures.
It automatically sets up a little filesystem for TRAL.
Moreover, you can decide which external software should be downloaded and installed without the need to read through their (sometimes complicated) installation procedures.


## Prerequisites

These installation scripts can be used on all Linux systems, but still require python3, pip, unzip and git.
Please make sure you use a python version >= Python3.5.
Per default TRAL will be installed within the site-packages of python and its external software in ``/usr/local/bin``. Then you may execute the scripts as root.
If you wish to change this path you can adjust 
``$INSTALLATION_PATH`` in ``configTRAL.cfg``

These setup scripts require ``bash``, ``python3``, ``pip``, ``unzip``, and ``git`` and run on any UNIX based operation system.
To compile some of the software (HMMER) ``gcc``, ``make`` and ``musl-dev`` is needed.
If you like to work in a virtual environment you need to install the respective software (e.g. virtenv)

## Getting Started

It is recommended to use TRAL within its own virtual environment for python3, e.g. with virtenv.

To install virtenv just use pip install:

```
pip3 install virtenv
```
Setup a virtenv and activate it:
```
virtualenv -p python3 env_tral
source env_tral/bin/activate
```
To get started with TRAL, simply clone the github directory https://github.com/acg-team/tral.git or download the directory easy_setup and then execute the bash scripts as described below. For example, create a folder called tral_repository in ``\home\user``:
```
mkdir tral_repository
cd tral_repository
git clone https://github.com/acg-team/tral.git
```

If you haven't git installed, you may need install git or to download the github repository manually (https://github.com/acg-team/tral.git).

## Installation of TRAL

You can adapt the default installing path within configTRAL_path.cfg. Please only change the variables FILES and ``$INSTALLATION_PATH``.

To install TRAL within you can simply run the script setupTRAL.sh with either the argument "pip" or "setup".
It is recommended to install TRAL within a virtual environment (e. g. virtualenv) which has to be activated before.

### Installation of TRAL with "setup"

Currently the version provided in "pypi" is not up-to-date, therefore it is recommended to install TRAL with python setup.py install "setup". 

```
cd ./tral/easy_setup
./setupTRAL.sh setup
```

Within the next release the tral version on pypi will be updated.

### Installation of external software for TRAL

After installing TRAL with setupTRAL.sh you can install external software. 
The script install_ext_software.sh automatically iterates through all installation scripts for each recommended external software of TRAL and ask you for each single software if you want to install it (press y/n). Please confirm to their respective licence.

```
sudo ./install_ext_software.sh
```

Otherwise, you can run an installation script for each external software individually which can be found within the directory ``setup_tral/install_ext_software``.
The software will be downloaded in the directory ``$FILES/tral_external_software``.
Currently, installation scripts for the following external software are available:

- alf
- hmmer
- hhrepid
- mafft
- phobos
- t-reks
- trf
- trust
- xstream

On the bottom of each installation script for the external software an uninstallation procedure can be found.
Either comment out all installation part and uncomment uninstallation part or run the commands directly in the commandline.
 
You may need to run these scripts as root (``$INSTALLATION_PATH`` is ``/usr/local/bin`` per default).

If you install these external software not in a user path like the default usr/local/bin, you have to adapt your config.ini file in .tral for each external software. The .tral directory will be copied to your ``$HOME`` directory automatically. TRAL only can find its configuration files if .tral is in your HOME directory. Otherwise provide a symlink (e.g. when installing on a cluster).

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
