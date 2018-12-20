.. _install_external:

Installation of external software
=================================

Here, we explain how to install external software packages, such as *de novo* tandem repeat
detectors. Important: Once a software is installed, the path to either the binary or the
executable shell script needs to be indicated in the TRAL configuration file :ref:`config.ini <configure>`.

Sequence profile model generation
---------------------------------

At current, only support for HMMER is integrated `published in Bioinformatics (2008) <http://bioinformatics.oxfordjournals.org/content/24/6/807.long>`_ (`Installation instructions <http://hmmer.janelia.org/>`__).
TRAL searches for HMMER's *hmmbuild* in the system path by default (modify in :ref:`config.ini <configure>`)::


    [hmm]
        hmmbuild = hmmbuild


If *hmmbuild* is not in your system path, set the absolute path::

    [hmm]
        hmmbuild = path/to/hmmbuild



.. _MAFFT:

Alignment of tandem repeat units
---------------------------------
Currently, MAFFT is the advised tool for (re-)alignment of the tandem repeat units to each other (`Installation instructions <http://mafft.cbrc.jp/alignment/software/>`__).
TRAL searches for MAFFT's *ginsi* in the system path by default (modify in :ref:`config.ini <configure>`)::

    [repeat]
        ginsi = ginsi


If *ginsi* is not in your system path, set the absolute path::

    [repeat]
        ginsi = path/to/ginsi


.. _install_denovo:

Currently integrated detectors
------------------------------

========    =======
Software    License
========    =======
HHrepID     CC-BY-NC-2.0
PHOBOS      Non-commercial
TRED        Non-commercial
T-REKS      Non-commercial?
TRF         Unlimited use
TRUST       ?
XSTREAM     Non-commercial
========    =======

HHrepID
*******

HHrepID is a profile self alignment *de novo* amino acid tandem repeat detection software
`published e.g. in PLOS Computational Biology (2011) <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195>`_ (`Installation instructions <http://toolkit.tuebingen.mpg.de/hhrepid>`__).
TRAL searches for the executable binary *hhrepid64* in the system path by default (modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            HHrepID = hhrepid_64

If the executable *hhrepid* is not in your system path, set the absolute path::

    [sequence]
        [[repeat_detector_path]]
            HHrepID = path/to/hhrepid

Also, you need to supply a null hmm file for using HHrepID. A dummy file is located in your home/.tral.
Supply the path to your null hmm file of choice::

    [sequence]
        [[repeat_detector_path]]
            HHrepID_dummyhmm = /path/to/home/.tral/data/dummyHMM.hmm


PHOBOS
******

`PHOBOS <http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm>`_  is an unpublished *k*-mer based *de novo* DNA tandem repeat detection software.
TRAL searches for the executable *phobos* in the system path by default (modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            PHOBOS = phobos


If *phobos* is not in your system path, set the absolute path::

    [sequence]
        [[repeat_detector_path]]
            PHOBOS = path/to/phobos


TRED
****

TRED is a sequence self alignment *de novo* amino acid tandem repeat detection software
`published in Bioinformatics (2007) <http://bioinformatics.oxfordjournals.org/content/23/2/e30.short>`_ (The software is available on request).
TRAL searches for the executable *tred* in the system path by default (modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            TRED = tred


If *phobos* is not in your system path, set the absolute path::

    [sequence]
        [[repeat_detector_path]]
            TRED = path/to/tred


T-REKS
******

T-REKS is a *k*-mer based *de novo* DNA/AA tandem repeat detection software
`published in Bioinformatics (2009) <http://bioinformatics.oxfordjournals.org/content/25/20/2632.short>`_ (`Installation instructions <http://bioinfo.montp.cnrs.fr/?r=t-reks>`__).
Create an executable text file T-REKS with the following content:
::

    #!/bin/sh
    java -jar /my/path/to/T-Reks.jar "$@"

If you place this text file in your systems path, TRAL finds it by default
(modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            T-REKS = T-REKS

If you did not place T-REKS in your system path or named it differently, set the absolute
path::

    [sequence]
        [[repeat_detector_path]]
            T-REKS = path/to/T-REKS


TRF
***

TRF is a *k*-mer based self alignment *de novo* DNA tandem repeat detection software
`published in Nucleic Acids Research (1999) <http://nar.oxfordjournals.org/content/27/2/573.full>`_ (`Installation instructions <http://tandem.bu.edu/trf/trf.html>`__).
TRAL searches for the executable *trf* in the system path by default (modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            TRF = trf


If *trf* is not in your system path, set the absolute path::

    [sequence]
        [[repeat_detector_path]]
            TRF = path/to/trf


TRUST
*****

TRUST is a sequence self alignment *de novo* amino acid tandem repeat detection software
`published in Bioinformatics (2004) <http://bioinformatics.oxfordjournals.org/content/20/suppl_1/i311.short>`_ (`Installation instructions <http://www.ibi.vu.nl/programs/trustwww/>`__).

Create an executable text file TRUST with the following content (you can amend the java
memory consumption restrictions)::

    #!/bin/sh
    java -Xmx30G -cp /my/path/to/TRUST/1.0.0/Align nl.vu.cs.align.SelfSimilarity "$@"

If you place this text file in your systems path, TRAL finds it by default
(modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            TRUST = TRUST

If you did not place TRUST in your system path or named it differently, set the absolute
path::

    [sequence]
        [[repeat_detector_path]]
            TRUST = path/to/TRUST

Also, you need to supply a substitution matrix for using TRUST (it ships with several substitution matrices).
Supply the path of your favourite substitution matrix::

    [sequence]
        [[repeat_detector_path]]
            TRUST_substitutionmatrix = /path/to/TRUST/Align/BLOSUM50


.. _XSTREAM:

XSTREAM
*******

XSTREAM is a *k*-mer based *de novo* DNA/AA tandem repeat detection software
`published in BMC Bioinformatics (2007) <http://www.biomedcentral.com/1471-2105/8/382/>`_ (`Installation instructions <http://jimcooperlab.mcdb.ucsb.edu/xstream/download.jsp>`__).

Create an executable text file XSTREAM with the following content:
::

    #!/bin/sh
    java -jar /my/path/to/xstream.jar "$@"

If you place this text file in your systems path, TRAL finds it by default
(modify in :ref:`config.ini <configure>`)::

    [sequence]
        [[repeat_detector_path]]
            XSTREAM = XSTREAM

If you did not place XSTREAM in your system path or named it differently, set the absolute
path::

    [sequence]
        [[repeat_detector_path]]
            XSTREAM = path/to/XSTREAM


Not yet integrated software
---------------------------

There is a large number of tandem repeat detection software for which TRAL does not provide
parsers. However, theses parsers are easily manually added to :ref:`sequence.repeat_detection_io <code_docs>`.
Please file an issue on the `tracker <https://github.com/elkeschaper/tandemrepeats/issues>`_.
