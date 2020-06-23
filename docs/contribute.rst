.. _contribute:

Technical hints for contributors
=================================


How to contribute
-----------------


1. Fork the repository
^^^^^^^^^^^^^^^^^^^^^^

`Forking <https://help.github.com/articles/fork-a-repo/>`_ a repository on github means creating a clone of a repository on github. Simply
click on "Fork" in the TRAL repos `TRAL repository <https://github.com/acg-team/tral/>`_
once you have a Github account.


2. Clone the repository from your fork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the command line, create a clone of your fork:

::

    git clone https://github.com/<YOUR_USER_NAME>/tral


3. Create your feature branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use git-flow, start a new feature:
::

    git-flow start feature <YOUR_FEATURE>


Otherwise, create a new branch as follows:
::

    git checkout -b feature/<YOUR_FEATURE>

Next, add the necessary changes and commit them:
::

    git add <CHANGED_FILE>
    git commit -m 'CHANGED_FILE: DESCRIPTION OF YOUR CHANGES’


4. Push your feature branch to your github fork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    git push origin feature/<FEATURE_NAME>


5. Create a pull request.
^^^^^^^^^^^^^^^^^^^^^^^^^

Create the `pull request <https://help.github.com/articles/using-pull-requests/>`_ online on github.
For this, go to your github page and click on "Pull Request".
::

    https://github.com/<YOUR_USER_NAME>/tral


Do a pull request on the develop branch of  acg-team/tral.

::

    base fork: acg-team/tral
    base: develop





How to contribute to the homepage
---------------------------------

Check out the current version of the TRAL homepage as follows:

::

    git clone --single-branch -b gh-pages https://<your_git_name>@github.com/acg-team/tral.git


How to test a package on pypitest
---------------------------------

Install the package from pypitest as follows:

::

    pip install -i https://testpypi.python.org/pypi



How to use Pypi
---------------


This documentation summarizes the steps to a new release of TRAL.

Version number updates

* README.md and README.rst
* setup.py
* docs/conf.py


Try testpypi

::
    cd ~
    swap .pypirc_testpypi .pypirc

    cd [TRAL directory]
    python setup.py sdist upload -r pypitest

    cd ~
    swap .pypirc_testpypi .pypirc

    cd [Installation test directory]
    pip install -i https://testpypi.python.org/pypi tral
