.. _background:

Background
==========

Printing tandem repeats
-----------------------
When printing a tandem repeat, you will see a lot of information. Here is what it means.
(You get this example from :ref:`the significance test tutorials <significance test>`) ::

    >>> print(lHIV_Sequence[4].dRepeat_list['denovo'].repeats[1])
    > begin:38 lD:2 n:6 pValue:0.0 divergence:0.46649169922545164 type:phylo_gap01
    RRN
    RR-
    RRW
    RA-
    RQ-
    RQI


The tandem repeat is displayed as an alignment of tandem repeat units, similar to
`multiple sequence alignments <http://en.wikipedia.org/wiki/Multiple_sequence_alignment>`_::

    RRN
    RR-
    RRW
    RA-
    RQ-
    RQI


Within the sequence, this tandem repeat would like as follows::

    RRNRRRRWRARQRQI


These are the other characteristics of the tandem repeat, that might be shown (if available):

- **begin**: This is the index where the tandem repeat starts in the sequence
- **lD**: This is the length of the consensus tandem repeat unit, not including insertions.
- **n**: This is the number of tandem repeat units.
- **pValue**: What is the probability that this repeat has occurred by random chance? (This value
  is a good estimate. So even if pValue=0.0, of course it is possible that the sequence shows similarity
  by random chance, and not because the repeat units have evolved by duplications.)
- **divergence**: this is a model-based measure of the similarity of the repeat units.
  (mathematically, it is the maximum likelihood estimate of the branch length on the phylogeny
  connecting all tandem repeat units. In the above example, every site has mutated 0.47 times
  on average. A value of 0 would mean no mutations have occurred, and the sequence of the
  tandem repeat units is very conserved.)
- **type**: Which model was used to calculate the pValue and the divergence?


More tandem repeat characteristics
----------------------------------

You have access to more tandem repeat characteristics. The dir() command will provide you
with a list of all attributes and functions connected to the object::

    >>> dir(lHIV_Sequence[4].dRepeat_list['denovo'].repeats[1])
    ['TRD', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'begin', 'calc_index_msa', 'calc_nD', 'calculate_pValues', 'calculate_scores', 'create', 'dDivergence', 'dPValue', 'dScore', 'deleteInsertionColumns', 'deletions', 'divergence', 'gapStructure', 'gap_structure_HMM', 'gaps', 'insertions', 'l', 'lD', 'msa', 'msaD', 'msaT', 'msaTD', 'msaTDN', 'msaTD_standard_aa', 'msa_original', 'msa_standard_aa', 'n', 'nD', 'nGap', 'pValue', 'save_original_msa', 'score', 'sequence_length', 'sequence_type', 'text', 'textD', 'textD_standard_aa', 'totD', 'write']


For a more explanation, check out the :ref:`Class documentation <code_docs>`, or
:ref:`contact us <contributors>`.

Writing a tandem repeat to .csv
--------------------------------

When you write a tandem repeat to .csv with TRAL, the result may look as follows::

    msa_original	lD	pValue	nD	sequence_length	begin
    GDII,GDIR	4	None	2.0	8	316
    FLG,FLG	3	None	2.0	6	507


Additional to the tandem repeat characteristics explained above, here you can find:

- **msa_original**: The tandem repeat unit alignment, with units separated by commata.
- **sequence_length**: The number of characters covered by the tandem repeat region.



