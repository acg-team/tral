.. _configure:

Configuration
=============

Here you learn how to configure TRAL after installation. All configuration files are located
in your home directory:
::

    ~/.tral/


configs.ini
------------

In this file, you need to define all paths and calculation defaults for using TRAL.


*De novo* repeat detection
**************************

Here, define the paths to all de novo detection algorithms you wish to use. The makeup of the
scripts or binaries needed is detailed in :ref:`Installation of external software <install_external>`.
Only paths for algorithms defined under [[repeat_detection]] are used, all other paths are
ignored.

::

    [sequence]
        [[repeat_detection]]
            AA = *List of default de novo detectors on amino acid data*
            DNA = *List of default de novo detectors on nucleic data*
        [[repeat_detector_path]]
            PHOBOS = path/to/PHOBOS
            HHrepID = path/to/HHREPID
            HHrepID_dummyhmm = path/to/HHREPID_HMM_DB
            TRUST = path/to/TRUST
            TRUST_substitutionmatrix = path/to/substitution_matrix/e.g.BLOSUM50 #Recommended: Use matrices packages with TRUST
            XSTREAM = path/to/XSTREAM
            [...]


Build Hmmer models
******************

::

    [hmm]
        hmmbuild = path/to/hmmbuild
        lDMax = 50  # maximum size of HMM for which the Viterbi algorithm is performed


Parameters of the default filter
********************************

::

    [filter]
        [[basic]]
            tag = basic_filter
            [[[dict]]]
                [[[[pValue]]]]
                    func_name = pValue
                    score = phylo_gap01  # score
                    threshold = 0.1  # p-Value cut-off
                [[[[nD]]]]
                    func_name = attribute
                    attribute = nD  # attribute of tandem repeat, e.g. the repeat unit length lD
                    type = min
                    threshold = 1.9  # p-Value cut-off



Default behaviour when new Repeat instances are created
*******************************************************

::

    [repeat]
        scoreslist = phylo_gap01,  # score
        calc_score = False  # Is the score calculated?
        calc_pValue = False # Is the pValue calculated?
        precision = 10


Default parameters for the model of repeat evolution used for statistical significance calculation
**************************************************************************************************

::

    [repeat_score]
        evolutionary_model = lg  # All models need to be available in /path/to/TandemRepeats/tandemrepeats/data/paml/
        [[indel]]
            indelRatePerSite = 0.01  # What magnitude is the indel rate compared to the substitution rate?
            ignore_gaps = True  # Shall trailing gaps be ignored / not penalised?
            gaps = row_wise
            zipf = 1.821
        [[optimisation]]
            start_min = 0.5
            start_max = 1.5
            nIteration = 14
        [[K80]]
            kappa = 2.59
        [[TN93]]
            alpha_1 = 0.3
            alpha_2 = 0.4
            beta = 0.7
    [[score_calibration]]
        scoreslist=phylo_gap01,
        save_calibration = False
        precision = 10


Default behaviour when new Repeat list instances are created
************************************************************

::

    [repeat_list]
        msa_original = True  # Is msa_original calculated?
        lD = True  # Is lD calculated?
        nD = True  # Is nD calculated?
        sequence_length = True  # Is the sequence_length calculated?
        pValue = phylo_gap01
        begin = True  # Is the position in sequence calculated?


logging.ini
-----------

In this file, you can define the level of debugging per module (DEBUG, INFO, WARNING), and
the format of the debugging message. Defaults to WARNING. The path to the file needs to be
defined as

::

    import logging
    import logging.config
    logging.config.fileConfig("path/to/your/home/.tral/logging.ini")



p-Value distribution files
--------------------------

In order to calculate the p-Value of tandem repeat scores, available p-Value distributions
need to be downloaded and placed in *./tral/data/pValue*:
::

    cd ~/.tral/data/pValue
    svn checkout https://github.com/elkeschaper/tral/trunk/tral/data/pValue .





