.. _install_containerized:

Using TRAL in a container or VM
===============================

Vagrant
-------

(`Vagrant <https://www.vagrantup.com/>`__) is tool for building virtual development
environments. TRAL requires a complex ecosystem to use all of its features,
so vagrant provides a reproducable environment to test installation.

Installing TRAL within a vagrant box is simple. We assume you have
(`installed vagrant <https://www.vagrantup.com/intro/getting-started/index.html>`__)
and a provider (e.g. VirtualBox). Next, check out the TRAL source code
and run
::

    vagrant up

This will download all dependencies and install all external software using
the :ref:`easy setup <easy_setup>` scripts.

You can then ssh into the box and work with TRAL.
::

    vagrant ssh
    $ sudo su
    # python3 -c 'import tral; print(tral.__version__)'

Note that the tral data files are downloaded to ``/root/.tral`` within the box,
so running (within vagrant) as root is recommended. The working directory
on the host is mounted as ``/vagrant`` within the box for easy data transfer.

Note that external software is licensed independently of TRAL. You should check
the licenses for all software. Several tools are restricted to academic use.
