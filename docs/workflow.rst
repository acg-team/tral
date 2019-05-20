.. _workflow:

Use GC3Pie to annotate your large sequence dataset.
===================================================

This tutorial will work you through an example of annotating tandem repeats on a very
large sequence data set: We will annotate tandem repeats with sequence profile models,
and with *de novo* detection algorithms, perform significance testing and overlap
filtering all in one.

For very large sequence sets (> 6000 sequences ~ 10h runtime), you may wish to perform the
annotation on several computing  nodes in parallel. Here, we provide a system to automise
the distribution and collection of tandem repeat annotation jobs with
`GC3Pie <https://code.google.com/p/gc3pie/>`_.

Everything is explained with a `toy workflow located in the TRAL source`_. You can print
the path to the local copy included in TRAL with this Python code::

    import os
    from tral.paths import PACKAGE_DIRECTORY
    print(PACKAGE_DIRECTORY)

To follow this tutorial, set a shell variable $MYTRAL to the package variable PACKAGE_DIRECTORY. In bash and sh, this can be done with::

    export MYTRAL=/path/to/your/package/directory

For csh or tcsh, use::

    setenv MYTRAL /path/to/your/package/directory

.. _`toy workflow located in the TRAL source`: https://github.com/acg-team/tral/tree/develop/tral/examples/workflow



Configure TRAL.
---------------

We start by :ref:`adapting the TRAL configuration files <configure>` to choose for
example which external tandem repeat detection algorithms you would like to run on the
sequence data. Don't forget to add a comma behind the last entry in this list, otherwise the
file format is violated.

Now, let's move ahead to data acquisition:


Prepare your data.
------------------
To annotate sequences with tandem repeats, we already acquired the following data and stored it in *$MYTRAL/examples/workflow* :

Sequence data
^^^^^^^^^^^^^^^^

We downloaded all `UniProtKB/Swiss-Prot human PRDM genes
<http://www.uniprot.org/uniprot/?query=gene%3Aprdm+AND+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score>`_
in fasta format...

::

    cd $MYTRAL/examples/workflow
    less uniprot_PRDM.fasta


... and split the sequence file into subsets on which the tandem repeat annotation can be
performed in parallel.

::

    ls split_sequence_data/*


On a sidenote, `fastasplitn <ftp://saf.bio.caltech.edu/pub/software/molbio/fastasplitn.c>`_
is a useful tool for splitting large data files.


Sequence profile model to sequence data mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We downloaded the mapping to PFAM annotations for these sequence also from www.uniprot.org.

::

    less uniprot_PRDM_annotation.tsv


Sequence profile models
^^^^^^^^^^^^^^^^^^^^^^^^

We downloaded all needed PFAM sequence profile models, for example of the `Zinc finger
<http://pfam.xfam.org/family/PF00096/hmm>`_.

::

    less PRDM_PFAM_models.hmm


pyfaidx
^^^^^^^^
The toy example requires the Python package `pyfaidx <https://github.com/mdshw5/pyfaidx>`_
for .fasta file indexing.

::

    pip install pyfaidx


Test the file preprocessing.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To allow fast access to the data for the tandem repeat annotation steps, it is transformed
to Python pickle objects:

::

    mkdir $MYTRAL/examples/workflow/hmm
    python3 tandem_repeat_annotation_scripts.py file_preparation --hmm_annotation_raw $MYTRAL/examples/workflow/uniprot_PRDM_annotation.tsv --hmm_annotation $MYTRAL/examples/workflow/uniprot_PRDM_annotation.pickle --hmm_raw $MYTRAL/examples/workflow/PRDM_PFAM_models.hmm --hmm $MYTRAL/examples/workflow/hmm


On the small toy example, this file preparation step should run fast.



Perform a test run.
-------------------
With the following command, annotation is performed on *uniprot_PRDM_1.fasta*:

::

    python3 tandem_repeat_annotation_scripts.py workflow -i $MYTRAL/examples/workflow/split_sequence_data/uniprot_PRDM_1.fasta -o $MYTRAL/examples/workflow/results/uniprot_PRDM_1.pickle -os $MYTRAL/examples/workflow/results/uniprot_PRDM_1.tsv -f tsv -t 600  --hmm_annotation $MYTRAL/examples/workflow/uniprot_PRDM_annotation.pickle --hmm $MYTRAL/examples/workflow/hmm


If this runs fine, you should see annotation results in:
::

    less results/uniprot_PRDM_1.tsv


Now, we can move ahead to automated distributed annotation with GC3PIE.


Install GC3PIE.
---------------

You can follow the official
`GC3Pie installation instructions <http://gc3pie.readthedocs.org/en/latest/users/install.html>`_.
Upon installation, GC3Pie needs to be `locally configured <http://gc3pie.readthedocs.org/en/latest/users/configuration.html>`_.


Usage
-----

Adapt this command to run the tandem repeat annotation workflow (`more information <http://gc3pie.readthedocs.org/en/latest/users/gc3apps/intro.html>`_)::


    $ ./tandem_repeat_annotation_workflow.py -w 60 minutes -r <host> -J 500 -u sqlite:////path/to/<session_name>.db -s <session_name> -C 2 -vvvv -conf $MYTRAL/examples/workflow/tandem_repeat_annotation_workflow.ini


Control the workflow with `GC3Utils <http://gc3pie.readthedocs.org/en/latest/users/gc3utils.html>`_.

::

    $ ./tandemrepeatannotationworkflow -s <session_name>
    $ gselect
    $ gresub
    $ gstat
    $ ...


