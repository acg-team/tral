.. _overlap_filtering:

Perform overlap filtering of redundant tandem repeat annotations.
=================================================================

Here you learn how to detect tandem repeat annotations from multiple sources for different
types of overlap, and filter them accordingly.

Requirements for this tutorial:

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



Filter overlapping tandem repeats.
----------------------------------------------------

Tandem repeats can be clustered according to different types of overlap:
::

    overlap_type = "shared_char"
    for iSequence in lHIV_Sequence:
        iSequence.dRepeat_list['denovo'].cluster(overlap_type)


In the first sequences, no repeats share any chars, however in the second sequence, three tandem repeats
overlap:
::

    >>> print(lHIV_Sequence[0].dRepeat_list['denovo'].dCluster[overlap_type])
    [{4}, {3}, {2}, {1}, {0}]
    >>> print(lHIV_Sequence[1].dRepeat_list['denovo'].dCluster[overlap_type])
    [{1, 2, 4}, {3}, {0}]


Tandem repeats can also directly be filtered of overlapping tandem repeats. Here, we need
to choose to retain one of the overlapping tandem repeats. For example the following usage
will first retain the tandem repeats of lowest p-Value according to the *score*, and in case
there still are draws, retain the tandem repeat of lowest divergence:
::

    overlap_type = "common_ancestry"
    score = "phylo_gap001"
    for iSequence in lHIV_Sequence:
        repeat_list_filtered = iSequence.dRepeat_list['denovo'].filter("none_overlapping", (overlap_type, None), [("pValue", score), ("divergence", score)])
        iSequence.set_repeat_list(repeat_list_filtered, "denovo_non_overlapping")


The resulting *result_list* now only contains tandem repeats that do not overlap according
the *common ancestry* overlap:
::

    >>> len([iR for iS in lHIV_Sequence for iR in iS.dRepeat_list["denovo"].repeats])
    22
    >>> len([iR for iS in lHIV_Sequence for iR in iS.dRepeat_list["denovo_non_overlapping"].repeats])
    18



Available types of overlap/redundancy detection
-----------------------------------------------

There are currently two types of overlap implemented:

- "shared_char": Do two tandem repeats contain any two chars in common?
- "common ancestry": Do two tandem repeats have at least two chars in the same column of their tandem repeat unit alignments? This approach is more conservative than "shared char". This approach has been used in the :ref:`MBE and New Phytologist publications, 2014 <publications>`.


For clustering, overlap is assumed to be a transitive attribute. That is, if tandem repeats
A & B, as well as B & C overlap, all tandem repeats A, B & C are clustered, no matter
whether A and C do also overlap.