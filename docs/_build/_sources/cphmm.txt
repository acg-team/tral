.. _cphmm:

Annotate tandem repeats from sequence domain models.
====================================================

Here you learn how to annotate sequences with your favourite sequence profile model.
The sequence profile model is transformed into a circular profile HMM, as described in
a :ref:`recent MBE publication <publications>`. Next, the sequence
is searched for tandem repeats using the circular profile HMM.

Requirements for this tutorial:

- :ref:`Install TRAL <install>`. TRAL ships with the data needed for this tutorial.
- :ref:`Install Mafft/ginsi <MAFFT>` for tandem repeat unit alignment.


Read in your sequence profile model.
------------------------------------

::

    import os
    from tral.hmm import hmm
    from tral.paths import PACKAGE_DIRECTORY

    pfam_profile_HMM = os.path.join(PACKAGE_DIRECTORY,"examples", "data","Kelch_1.hmm")
    circular_profile_HMM_Kelch_1 = hmm.HMM.create(input_format = 'hmmer', file = pfam_profile_HMM)


See the basic characteristics of the circular profile HMM::

    >>> print(circular_profile_HMM_Kelch_1)
    cpHMM ID: PF01344.20
    cpHMM length: 47
    Most likely motif: ARSSAGVVVLDGKIYVIGGRDGDGNALNSVERYDPVTNTWEKLPSMP


Read in your sequences.
-----------------------

::

    from tral.sequence import sequence

    human_HCFC1_fasta = os.path.join(PACKAGE_DIRECTORY,"examples", "data", "P51610.fasta")
    human_HCFC1_sequence = sequence.Sequence.create(file = human_HCFC1_fasta, input_format = 'fasta')[0]




Annotate tandem repeats with the circular profile HMM.
------------------------------------------------------

::

    tandem_repeats = human_HCFC1_sequence.detect(lHMM = [circular_profile_HMM_Kelch_1])


The result is a tandem repeat :ref:`(interpretation) <background>`::

    >>> print(tandem_repeats.repeats[0])
    > begin:31 l_effective:53 n:5
    R--------PRHGHRAVAIKELIVVFGGGN----------EGIVD-----------------------------------------------------------ELHVYNTATNQWFI---PAVRGDIP-
    P--------GCAAYGFVCDGTRLLVFGGMV-----------------------------------EYG-------------------KYSN-------------DLYELQASRWEWKR-----LKAK---
    TPKNGPPPCPRLGHSFSLVGNKCYLFGGLANDSEDPKNNIPRYLNDLYILELRPGSGVVAWDIPITYGVLPPPRESHTAVVYTEKDNKKSKLVIYGGMSGCRLGDLWTLDIDTLTWNK---PSLSGVAPL
    ---------PRSLHSATTIGNKMYVFGGWV----------PLVMDDV-------------------------------KVATHEKEWKCTN-------------TLACLNLDTMAWETILMDTLEDNIP-
    R--------ARAGHCAVAINTRLYI---------------------------------------------------------------------------------------------------------

Output the detected tandem repeats.
-----------------------------------

Write a singe repeat_list to .tsv format::

    path_to_output_tsv_file = "outputfile.tsv"  # Choose your path and filename
    tandem_repeats.write(output_format = "tsv", file = path_to_output_tsv_file)


This is how the output looks like :ref:`(interpretation) <background>`::

    $ cat outputfile.tsv
    begin	msa_original	l_effective	n_effective	repeat_region_length	divergence	pvalue
    31	R--------PRHGHRAVAIKELIVVFGGGN----------EGIVD-----------------------------------------------------------ELHVYNTATNQWFI---PAVRGDIP-,P--------GCAAYGFVCDGTRLLVFGGMV-----------------------------------EYG-------------------KYSN-------------DLYELQASRWEWKR-----LKAK---,TPKNGPPPCPRLGHSFSLVGNKCYLFGGLANDSEDPKNNIPRYLNDLYILELRPGSGVVAWDIPITYGVLPPPRESHTAVVYTEKDNKKSKLVIYGGMSGCRLGDLWTLDIDTLTWNK---PSLSGVAPL,---------PRSLHSATTIGNKMYVFGGWV----------PLVMDDV-------------------------------KVATHEKEWKCTN-------------TLACLNLDTMAWETILMDTLEDNIP-,R--------ARAGHCAVAINTRLYI---------------------------------------------------------------------------------------------------------	53	4.056603773584905	306	None	None


Write a singe repeat_list to .pickle format::

    path_to_output_pickle_file = "outputfile.pickle"  # Choose your path and filename
    tandem_repeats.write(output_format = "pickle", file = path_to_output_pickle_file)


A repeat_list in pickle format can easily be read in again::

    from tral.repeat_list import repeat_list
    tandem_repeats = repeat_list.RepeatList.create(input_format = "pickle", file = path_to_output_pickle_file)