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

* README.md
* tral/__init__.py
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

Release Process
===============
1. Bump Version Number
----------------------
Change version in:

* README.md
* tral/__init__.py
* docs/conf.py

2. Test  
-----------------
Run all tests again.
::

    cd tral/tral
    pytest

Make sure, that `pytest` is installed in the virtual environment.

3. Rebuild docs
-----------------
The documentation is build with `sphinx`.
Check if all external links are correct, update and commit if necessary.

Checkout branch develop:
::

    git checkout develop

Move into './tral/docs' directory and build the docs:
::

    cd tral/docs
    make html

The new docs are stored in `./tral/docs/_build/html/`.
Copy them outside the git repository and checkout to the `gh-pages` branch.
::

    cp -r tral/docs/_build/html ~
    git checkout gh-pages

Move the fresh build docs to `gh-pages`:
::

    mv -r ~/html/* ./tral/docs

[CHECK IF THIS IS CORRECT]

Check the updated docs online (synchronizes automatically) at `acg-team.github.io <https://acg-team.github.io/tral/index.html>`_.

4. Build a python wheel
-----------------------
::

    git checkout develop  
    python setup.py sdist bdist_wheel

5. Create a git tag
-------------------
`-s` requires a `gpg key <https://git-scm.com/book/en/v2/Git-Tools-Signing-Your-Work>`_.
::

    git tag -s -a <VERSION>

6. Upload to PyPi test
----------------------
Check the wheel and git tag:
::

    twine check dist/*; git verify-tag <VERSION>

Upload to PyPi test platform first:
::

    twine upload --repository testpypi dist/*

7. Git Release
--------------
::

    git checkout develop
    git push; git push --tags

8. Upload to PyPi
-----------------
::
    
    twine upload dist/*

8. Build docker container
--------------------------
Create `Github.com Token <https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line>`_
and login your docker instance:

::
    
    cat TOKEN.txt | docker login https://docker.pkg.github.com -u USERNAME --password-stdin
    
We encountered issues if not all tral related docker container and images were stopped and removed.
Continue then with building a new docker image with the new release:

::
    
    sudo docker build -t docker.pkg.github.com/acg-team/tral/tral_docker:<VERSION> -t docker.pkg.github.com/acg-team/tral/tral_docker:latest -f tral/docker/Dockerfile --no-cache 
    sudo docker push docker.pkg.github.com/acg-team/tral/tral_docker:latest


10. Prepare repository for further development
----------------------------------------------
Increment version in `__init_.py` from i.e. `2.0.0` to `2.0.1.dev0`.
