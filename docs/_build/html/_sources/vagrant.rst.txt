.. _vagrant:

Running TRAL using vagrant
==========================

Building
--------

To build a vagrant box, first clone tral from the git repo. Next, run
::

    vagrant up
    vagrant ssh

This will The current working directory will be mounted at /vagrant within the box.
Tral should be run as root within the box (it uses /root/.tral as the data dir).
::

    sudo su
    cd /vagrant


Testing TRAL
------------

The vagrant box is useful for testing tral.
::

    vagrant ssh
    sudo su
    cd /vagrant
    pytest
