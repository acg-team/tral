.. _install_external:

Installation of external software
=================================

Here, we explain how to install external software packages, such as *de novo* tandem repeat
detectors. Important: Once a software is installed, the path to either the binary or the
executable shell script needs to be indicated in the TRAL configuration file :ref:`defaults.ini <configure>`.



Sequence profile model generation
---------------------------------

At current, only support for HMMER is integrated `published in Bioinformatics (2008) <http://bioinformatics.oxfordjournals.org/content/24/6/807.long>`_ (`Installation instructions <http://hmmer.janelia.org/>`__).

::

    hmmbuild = /my/path/to/hmmbuild


Alignment of tandem repeat units
---------------------------------
Currently, MAFFT is the advised tool for (re-)alignment of the tandem repeat units to each other (`Installation instructions <http://mafft.cbrc.jp/alignment/software/>`__).

::

    ginsi = /my/path/to/ginsi



Currently integrated detectors
------------------------------

HHrepID
*******

HHrepID is a profile self alignment *de novo* amino acid tandem repeat detection software
`published e.g. in PLOS Computational Biology (2011) <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195>`_ (`Installation instructions <http://toolkit.tuebingen.mpg.de/hhrepid>`__).

The executable is the binary *hmmbuild*:
::

    /my/path/to/hhrepid_64


PHOBOS
******

`PHOBOS <http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm>`_  is an unpublished *k*-mer based *de novo* DNA tandem repeat detection software.

The executable is the binary:
::

    /my/path/to/phobos


TRED
****

TRED is a sequence self alignment *de novo* amino acid tandem repeat detection software
`published in Bioinformatics (2007) <http://bioinformatics.oxfordjournals.org/content/23/2/e30.short>`_ (The software is available on request).

The executable is the binary:
::

    /my/path/to/tred



T-REKS
******

T-REKS is a *k*-mer based *de novo* DNA/AA tandem repeat detection software
`published in Bioinformatics (2009) <http://bioinformatics.oxfordjournals.org/content/25/20/2632.short>`_ (`Installation instructions <http://bioinfo.montp.cnrs.fr/?r=t-reks>`__).

Create an executable text file with the following content:
::

    #!/bin/sh
    java -jar /my/path/to/T-Reks.jar "$@"


TRF
***

TRF is a *k*-mer based self alignment *de novo* DNA tandem repeat detection software
`published in Nucleic Acids Research (1999) <http://nar.oxfordjournals.org/content/27/2/573.full>`_ (`Installation instructions <http://tandem.bu.edu/trf/trf.html>`__).

The executable is the binary:
::

    /my/path/to/trf


TRUST
*****

TRUST is a sequence self alignment *de novo* amino acid tandem repeat detection software
`published in Bioinformatics (2004) <http://bioinformatics.oxfordjournals.org/content/20/suppl_1/i311.short>`_ (`Installation instructions <http://www.ibi.vu.nl/programs/trustwww/>`__).

Create an executable text file with the following content (you can amend the java memory consumption restrictions):
::

    #!/bin/sh
    java -Xmx30G -cp /my/path/to/TRUST/1.0.0/Align nl.vu.cs.align.SelfSimilarity "$@"


.. _XSTREAM:

XSTREAM
*******

XSTREAM is a *k*-mer based *de novo* DNA/AA tandem repeat detection software
`published in BMC Bioinformatics (2007) <http://www.biomedcentral.com/1471-2105/8/382/>`_ (`Installation instructions <http://jimcooperlab.mcdb.ucsb.edu/xstream/download.jsp>`__).

Create an executable text file with the following content:
::

    #!/bin/sh
    java -jar /my/path/to/xstream.jar "$@"



Not yet integrated software
---------------------------

There is a large number of tandem repeat detection software for which TRAL does not provide
parsers. However, theses parsers are easily manually added to :ref:`sequence.repeat_detection_io <sequence>`.
Please file an issue on the `tracker <https://github.com/elkeschaper/tandemrepeats/issues>`_.