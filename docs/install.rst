.. _install:

Installation
============

Try:
::

    $ pip install tral


The TRAL configuration files are automatically added to your home directory:
::

    ~/.tral/


Install with sudo
-----------------

In case you install TRAL with sudo, the installation home might differ for the sudo account.
You can either supply the path to your home directory during installation ...::

    $ python setup.py install --home /path/to/your/home


or mv the TRAL configuration directory from your root dir to your home, e.g.::

    $ cd ~
    $ sudo mv /root/.tral ./
    $ sudo chown -R $(whoami) .tral
