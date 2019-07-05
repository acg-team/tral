.. _tral:

Tandem Repeat Annotation Library
================================

TRAL makes annotation of tandem repeats in amino acid and nucleic data simple. TRAL includes
modules for detecting tandem repeats with both de novo software and sequence profile HMMs;
statistical significance analysis of putative tandem repeats, and filtering of redundant predictions.


Getting started
===============

A working TRAL installation consists of several parts:

1. The TRAL python package
2. External tools and repeat detectors.
3. TRAL configuration files and supporting data.

TRAL itself requires Python3 and should run on all platforms. However, many
external tools require linux. The full TRAL pipeline has been tested on ubuntu.

The easiest way to get started on ubuntu is to use the
:ref:`easy setup system <easy_setup>`. This automatically installs TRAL and all
external tools.

Alternately, follow these steps to install each part:

1. :ref:`Install TRAL <install>`
2. :ref:`Install external software needed for TRAL <install_external>`
3. :ref:`Download p-value distribution files <pvaluefiles>`

Finally, :ref:`adapt the configuration files <configure>`.

If you have trouble installing TRAL on your system, consider using the
:ref:`containerized install <install_containerized>` .

Tutorials
=========

- :ref:`Run and parse <denovo>` de novo repeat detection software.
- Background about :ref:`tandem repeat characteristics in TRAL <background>`
- :ref:`Annotate tandem repeats<cphmm>` from sequence domain models.
- :ref:`Perform statistical significance test <significance_test>` of tandem repeats.
- :ref:`Perform overlap filtering <overlap_filtering>` of redundant tandem repeat annotations.
- :ref:`Use GC3Pie <workflow>` to annotate your large sequence dataset (toy example)
- :ref:`Find instances of a particular repeat<search_hmm>` in a sequence database.


Reference Guide
===============

- `Code base on GitHub <https://github.com/acg-team/tral>`_
- :ref:`Class documentation <code_docs>`
- `GC3Pie <https://code.google.com/p/gc3pie/>`_



How to cite us
===============

E Schaper, A Korsunsky, J Pecerska, A Messina, R Murri, H Stockinger, S Zoller, I Xenarios, and M Anisimova (2015). `TRAL: Tandem Repeat Annotation Library <http://bioinformatics.oxfordjournals.org/content/early/2015/05/17/bioinformatics.btv306.abstract>`_. *Bioinformatics*. DOI:  10.1093/bioinformatics/btv306

More on :ref:`contributors and the background <contributors>` of this project.


How to contribute
==================


- :ref:`Technical hints <contribute>` for contributors to TRAL.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* `search </search.html>`_

.. toctree::
   :hidden:

   install
   code_docs
   contribute
   contributors
   denovo
   background
   cphmm
   significance_test
   overlap_filtering
   workflow
   search_hmm

