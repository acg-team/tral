.. _denovo:

Run and parse de novo repeat detection software.
================================================

Here you learn how to run external repeat detection software on your favourite sequence
file and output the results.

Requirements for this tutorial:

- :ref:`Install TRAL <install>`.
- Install :ref:`XSTREAM <XSTREAM>`. If preferred, install one or more other :ref:`tandem repeat detectors<install_denovo>` instead.


Read in your sequences.
-----------------------

::

    import os
    from tral.sequence import sequence
    from tral.paths import PACKAGE_DIRECTORY

    fHIV_proteome = os.path.join(PACKAGE_DIRECTORY,"test","HIV-1_388796.faa")
    lHIV_Sequence = sequence.Sequence.create(file = fHIV_proteome, format = 'fasta')



Run the external repeat detector.
---------------------------------

Detect tandem repeats on a single sequence with :ref:`XSTREAM <XSTREAM>`. If you did not
install XSTREAM, you can use any of the other *de novo* detection algorithms, but will see
different results.
::

    tandem_repeats = lHIV_Sequence[0].detect(denovo = True, detection = {"lFinders": ["XSTREAM"]})


As an example, the first detected putative tandem repeat looks as follows::

    >>> print(tandem_repeats.repeats[0])
    > begin:316 lD:4 n:2
    GDII
    GDIR



Detect tandem repeats on all sequences with all *de novo* tandem repeat detection algorithms
defined in the :ref:`configuration file <configure>`::

    for iSequence in lHIV_Sequence:
        iTandem_repeats = iSequence.detect(denovo = True)
        iSequence.set_repeat_list(iTandem_repeats, "denovo")


As an example, the first detected putative tandem repeat looks as follows:
::

    >>> print(lHIV_Sequence[0].dRepeat_list['denovo'].repeats[0])
    > begin:31 lD:53 n:5
    R--------PRHGHRAVAIKELIVVFGGGN----------EGIVD-----------------------------------------------------------ELHVYNTATNQWFI---PAVRGDIP-
    P--------GCAAYGFVCDGTRLLVFGGMV-----------------------------------EYG-------------------KYSN-------------DLYELQASRWEWKR-----LKAK---
    TPKNGPPPCPRLGHSFSLVGNKCYLFGGLANDSEDPKNNIPRYLNDLYILELRPGSGVVAWDIPITYGVLPPPRESHTAVVYTEKDNKKSKLVIYGGMSGCRLGDLWTLDIDTLTWNK---PSLSGVAPL
    ---------PRSLHSATTIGNKMYVFGGWV----------PLVMDDV-------------------------------KVATHEKEWKCTN-------------TLACLNLDTMAWETILMDTLEDNIP-
    R--------ARAGHCAVAINTRLYI---------------------------------------------------------------------------------------------------------


Output the detected tandem repeats.
-----------------------------------

Write a singe repeat_list to .tsv format::

    path_to_output_tsv_file = "/my/path/to/the/outputfile.tsv"
    tandem_repeats.write(format = "tsv", file = path_to_output_tsv_file)


The created .tsv looks as follows::

    $ cat /my/path/to/the/outputfile.tsv
    msa_original	lD	pValue	nD	sequence_length	begin
    GDII,GDIR	4	None	2.0	8	316
    FLG,FLG	3	None	2.0	6	507


Write a singe repeat_list to .pickle format::

    path_to_output_pickle_file = "/my/path/to/the/outputfile.pickle"
    tandem_repeats.write(format = "pickle", file = path_to_output_pickle_file)


A repeat_list in pickle format can easily be read in again::

    from tral.repeat_list import repeat_list
    tandem_repeats = repeat_list.Repeat_list.create(format = "pickle", file = path_to_output_pickle_file)


Save multiple sequence together with tandem repeat annotations::

    import pickle
    path_to_output_pickle_file = "/my/path/to/the/outputfile.pickle"
    with open(path_to_output_pickle_file, 'wb') as fh:
        pickle.dump(lHIV_Sequence, fh)
