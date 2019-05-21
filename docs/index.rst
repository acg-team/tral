.. tral

Tandem Repeat Annotation Library
================================

TRAL makes annotation of tandem repeats in amino acid and nucleic data simple. TRAL includes
modules for detecting tandem repeats with both de novo software and sequence profile HMMs;
statistical significance analysis of putative tandem repeats, and filtering of redundant predictions.


Getting started
===============

To get started with TRAL you need to follow several steps.

1. Installation of all relevant software and files:

	* 1(a-c) Simplified setup with the :ref:`easy setup system for TRAL <easy_setup>`.

   Otherwise download and install everything manually:

	* 1(a) :ref:`Install TRAL <install>`
	* 1(b) :ref:`Install external software needed for TRAL <install_external>`
	* 1(c) :ref:`Download p-value distribution files <pvaluefiles>`

2. :ref:`Adapt configuration files <configure>`


Tutorials
=========

- :ref:`Run and parse <denovo>` de novo repeat detection software.
- Background about :ref:`tandem repeat characteristics in TRAL <background>`
- :ref:`Annotate tandem repeats<cphmm>` from sequence domain models.
- :ref:`Perform statistical significance test <significance_test>` of tandem repeats.
- :ref:`Perform overlap filtering <overlap_filtering>` of redundant tandem repeat annotations.
- :ref:`Use GC3Pie <workflow>` to annotate your large sequence dataset (toy example)


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
* :ref:`search`
