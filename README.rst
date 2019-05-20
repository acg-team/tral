Tandem Repeat Annotation Library
================================

TRAL is a highly modularized, flexible sequence tandem repeats
annotation Python2/3 library.

-  Large scale annotation of tandem repeats with *de novo* detectors,
   and sequence profile models
-  Statistical significance testing, overlap detection, and filtering of
   annotations
-  Refinement of tandem repeat annotations with circular profile hidden
   Markov models
-  User-defined output formats

The source code is `documented on GitHub
IO <http://acg-team.github.io/tral/>`__.

Version
~~~~~~~

0.3.5

Installation
~~~~~~~~~~~~

TRAL is available on `Pypi <https://pypi.python.org/pypi>`__ and can be
installed with `pip <https://pip.pypa.io/en/latest/>`__ for Python>=3.2:

.. code:: sh

    $ pip install tral

See also more extensive `Installation
instructions <http://acg-team.github.io/tral/install.html#install>`__.

License
~~~~~~~

GPL-2.0

Dependencies
~~~~~~~~~~~~

Some of TRAL's functions depend on external software (`Installation instructions for dependencies <http://acg-team.github.io/tral/install_external.html#install-external>`__). This includes creation of sequence profile hidden Markov models, alignment of tandem repeat units, and *de novo* repeat detection.

