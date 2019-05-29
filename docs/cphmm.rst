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
- :ref:`Install Castor`. If you wish to realign your repeat MSA you need to have proPIP installed beforehand.

.. todo:: write section to install castor with aligner proPIP


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

Annotate tandem repeats with the circular profile HMM and realign with proPIP.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Alternatively, you can directly include tandem repeat MSA realignment with the indel-aware proPIP algorithm::

	tandem_repeats = human_HCFC1_sequence.detect(lHMM = [circular_profile_HMM_Kelch_1], realignment='proPIP')

::

	>>> print(tandem_repeats.repeats[0])
	> begin:31 l_effective:48 n:5
	-RP-R-----------HGHRA-V---A--------I-----K--------ELIVVF----G-G----------GN-----E-G---I--V--D--ELH-V-YN-T-A-T--N--Q--W---F---IPAV---R-GD---I-P-
	--PGC-----A---A-YGF---V---C--D-G---T-----R---------LL-VF----G-G-------M---V-----EYG--KY--S--N--DL----YELQ-A-S-----R--WE--W-KRLKA----K----------
	T-P-KNGPPPCPR-LGHSF-SLVGNKCYL-FGGLANDSEDPKNNIPRYLNDLY-ILELRPGSGVVAWDIPITYGVLPPPRE-SHTAVVYTEKDNKKSKLVIYG-GMSGCRLG-D-L-W-T-LD--IDTLTWNKP-SLSGVAPL
	--P-R-----S---L-HS--A-T---T-I--G---N-----K---------MY-VF----G-G---W-VPL---VM----D-D-VKVA-T-HE-KEWK-C-TN-TLA-C-LNLDTMAWETIL---MDTL---E--D-N-I-P-
	----R-----A--RA--GH-------C--------A------------------V-------A-----------I-----------------N---------------T-----R--L---Y---I-----------------


If you wish to use a gamma distribution for indels, choose the option rate_distribution = 'gamma' or modify the rate distribution in :ref:`config.ini <configure>`.

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
