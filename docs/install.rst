.. _install:

Installation
============


With pip (or pip3) configured for Python3, try::

    $ pip install tral


The TRAL configuration files are automatically added to your home directory:
::

    ~/.tral/


Install with sudo
-----------------

In case you install TRAL with sudo, the installation home might differ for the sudo account.
You can either supply the path to your home directory during installation ...::

    $ python setup.py install --tral-home /path/to/your/home


or mv the TRAL configuration directory from your root dir to your home, e.g.::

    $ cd ~
    $ sudo mv /root/.tral ./
    $ sudo chown -R $(whoami) .tral


Development
-----------

To compile locally, run

    $ pip install -r requirements_dev.py
    $ python setup.py install

TRAL uses tox and pytest for testing. The following commands relate to testing and code style:

    $ flake8 tral  # check code style
    $ python setup.py test  # Run tests
    $ tox  # Run comprehensive test suite on multiple python versions
