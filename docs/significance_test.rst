
.. _significance_test:

Perform statistical significance test of tandem repeats.
========================================================

Here you learn how to perform statistical significance tests on your candidate set of tandem repeats.
TRAL currently provides a variety of model-based and none-model based tests. They are
detailed in a :ref:`NAR publication (2012) <publications>`.

The requirements for this tutorial are:

- :ref:`Install TRAL <install>`.
- :ref:`Download p-Value distribution files <pValuefiles>`.


Read in tandem repeat annotations.
----------------------------------

::

    import os, pickle
    from tral import sequence
    from tral.paths import PACKAGE_DIRECTORY

    fRepeat_Pickle = os.path.join(PACKAGE_DIRECTORY,"test","HIV-1_388796.pickle")

    with open(fRepeat_Pickle, 'rb') as fh:
        lHIV_Sequence = pickle.load(fh)

    print(lHIV_Sequence)



Perform a statistical significance test.
----------------------------------------

Calculate the statistical score on the 'denovo' *repeat_list*, and the corresponding
p-Value with parameters defined in the :ref:`configuration file <configure>`:
::

    for iSequence in lHIV_Sequence:
        if iSequence.dRepeat_list['denovo']:
            for iTandemRepeat in iSequence.dRepeat_list['denovo'].repeats:
                iTandemRepeat.calculate_pValues()

For example, the following putative tandem repeat is found to be non-significant with the used model
of tandem repeat evolution:
::

    >>> print(lHIV_Sequence[2].dRepeat_list['denovo'].repeats[0])
    >begin:94 lD:7 n:2 pValue:1.0 divergence:1.656005859375
    EKGGLEGLIYSKKRQE
    ---ILDLWVY------


Whereas this putative tandem repeat is considered significant:
::

    >>> print(lHIV_Sequence[4].dRepeat_list['denovo'].repeats[1])
    >begin:38 lD:2 n:6 pValue:0.0 divergence:0.46649169922545164
    RRN
    RR-
    RRW
    RA-
    RQ-
    RQI


Besides, parameters can also be directly supplied to *calculate_pValues*:
::

    for iSequence in lHIV_Sequence:
        if iSequence.dRepeat_list['denovo']:
            for iTandemRepeat in iSequence.dRepeat_list['denovo'].repeats:
                iTandemRepeat.calculate_pValues(scoreslist=['phylo_gap001'])




Filter tandem repeat below a significance threshold.
----------------------------------------------------

::

    for iSequence in lHIV_Sequence:
        repeat_list = iSequence.dRepeat_list["denovo"]
        if repeat_list:
            repeat_list_filtered = repeat_list.filter(func_name = "pValue", score = "phylo_gap01", threshold = 0.05)
            iSequence.set_repeat_list(repeat_list_filtered, "denovo_filtered")

The resulting *result_list* now only contains tandem repeats with a p-Value below
0.05, and is shorter than the original tandem repeat list:

::

    >>> len([iR for iS in lHIV_Sequence for iR in iS.dRepeat_list["denovo"].repeats])
    22
    >>> len([iR for iS in lHIV_Sequence for iR in iS.dRepeat_list["denovo_filtered"].repeats])
    17



Available significance tests
----------------------------

The model based scores are generally described in our :ref:`NAR publication, 2012 <publications>`.
Briefly, the sequences evolution is assumed to be described by the LG matrix.

Further, the models can be distinguished by how gaps are handled:

- 'phylo': Gaps are treated as characters.

Gaps can also assumed to be  inserted with exponentially distributed waiting times and Zipfian distributed indel lengths:

- 'phylo_gap01', 'phylo_gap01_ignore_trailing_gaps', 'phylo_gap01_ignore_coherent_deletions',
  'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions': Gaps are modelled, and assumed to mutation at at 0.1 times lower rate compared to substitutions.
- 'phylo_gap001', 'phylo_gap001_ignore_trailing_gaps', 'phylo_gap001_ignore_coherent_deletions',
  'phylo_gap001_ignore_trailing_gaps_and_coherent_deletions': Gaps are modelled, and assumed to mutation  at at 0.01 times lower rate compared to substitutions.
-  *ignore_trailing_gaps* signifies that gaps before the first tandem repeat unit and after the last tandem repeat unit are not penalised.
-  *coherent_deletions* signifies that gaps of the same length and position within the alignment of tandem repeat units are penalised only once.

Among the three *ad hoc* scores are 'entropy', 'parsimony', and 'pSim' scores
(see our :ref:`NAR publication, 2012 <publications>`).

