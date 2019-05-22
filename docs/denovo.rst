.. _denovo:

Run and parse de novo repeat detection software.
================================================

Here you learn how to run external repeat detection software on your favourite sequence
file and output the results.

Requirements for this tutorial:

- :ref:`Install TRAL <install>`. TRAL ships with the data needed for this tutorial.
- :ref:`Install XSTREAM <XSTREAM>`. You can also install one or more other :ref:`tandem repeat detectors<install_denovo>` instead.


Read in your sequences.
-----------------------

::

    import os
    from tral.sequence import sequence
    from tral.paths import PACKAGE_DIRECTORY

    proteome_HIV = os.path.join(PACKAGE_DIRECTORY, "examples", "data", "HIV-1_388796.faa")
    sequences_HIV = sequence.Sequence.create(file = proteome_HIV, input_format = 'fasta')



Run the external repeat detector.
---------------------------------

Detect tandem repeats on a single sequence with :ref:`XSTREAM <XSTREAM>`. If you did not
install XSTREAM, you can use any of the other *de novo* detection algorithms, but will see
different results.
::

    tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["XSTREAM"]})


As an example, the first detected putative tandem repeat looks as follows :ref:`(interpretation) <background>`::

    >>> print(tandem_repeats.repeats[0])
    > begin:316 l_effective:4 n:2
    GDII
    GDIR



Detect tandem repeats on all sequences with all *de novo* tandem repeat detection algorithms
defined in the :ref:`configuration file <configure>`. Here "denovo" is used as a tag, you can choose your own tag.::

    for iSequence in sequences_HIV:
        iTandem_repeats = iSequence.detect(denovo = True)
        iSequence.set_repeatlist(iTandem_repeats, "denovo")


Different different algorithms usually detect different tandem repeats. This the the
absolute number of detections in the HIV proteome for a couple of algorithms::

    >>> len([i for j in sequences_HIV for i in j.get_repeatlist('denovo').repeats if i.TRD == "HHrepID"])
    9
    >>> len([i for j in sequences_HIV for i in j.get_repeatlist('denovo').repeats if i.TRD == "T-REKS"])
    4
    >>> len([i for j in sequences_HIV for i in j.get_repeatlist('denovo').repeats if i.TRD == "TRUST"])
    6
    >>> len([i for j in sequences_HIV for i in j.get_repeatlist('denovo').repeats if i.TRD == "XSTREAM"])
    9


As an example, T-REKS detects the following repeat in the second HIV sequence :ref:`(interpretation) <background>`::

    >>> print([i for i in sequences_HIV[1].get_repeatlist('denovo').repeats if i.TRD == "T-REKS"][0])
    > begin:449 l_effective:10 n:2
    RPEPTAPP-ESL
    RPEPTAPPPES-

Output the detected tandem repeats.
-----------------------------------

Write a singe repeat_list to .tsv format::

    path_to_output_tsv_file = "outputfile.tsv" # Choose your path and filename
    tandem_repeats.write(output_format = "tsv", file = path_to_output_tsv_file)


The created .tsv looks as follows :ref:`(interpretation) <background>`::

    $ cat outputfile.tsv
    begin   msa_original    l_effective      nD      sequence_length divergence      pvalue
    316     GDII,GDIR       4       2.0     8       None    None
    507     FLG,FLG 3       2.0     6       None    None


Write a singe repeat_list to .pickle format::

    path_to_output_pickle_file = "outputfile.pickle"  # Choose your path and filename
    tandem_repeats.write(output_format = "pickle", file = path_to_output_pickle_file)


A repeat_list in pickle format can easily be read in again::

    from tral.repeat_list import repeat_list
    tandem_repeats = repeat_list.RepeatList.create(input_format = "pickle", file = path_to_output_pickle_file)


Save multiple sequence together with tandem repeat annotations::

    import pickle
    path_to_output_pickle_file = "outputfile.pickle" # Choose your path and filename
    with open(path_to_output_pickle_file, 'wb') as fh:
        pickle.dump(sequences_HIV, fh)
