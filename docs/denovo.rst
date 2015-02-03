.. _denovo:

Run and parse de novo repeat detection software.
================================================

Here you learn how to run external repeat detection software on your favourite sequence
file and output the results.


Read in your sequences.
-----------------------

::

    import os
    from tandemrepeats.sequence import sequence

    path_to_tandemrepeats = "/my/path/to/TandemRepeats/tandemrepeats"

    fHIV_proteome = os.path.join(path_to_tandemrepeats,"test","HIV-1_388796.faa")
    lHIV_Sequence = sequence.Sequence.create(file = fHIV_proteome, format = 'fasta')



Run the external repeat detector.
---------------------------------


Detect tandem repeats on a single sequence with :ref:`XSTREAM <XSTREAM>`:
::

    tandem_repeats = lHIV_Sequence[0].detect(denovo = True, detection = {"lFinders": ["XSTREAM"]})
    print(tandem_repeats.repeats)
    print(tandem_repeats.repeats[0])


Detect tandem repeats on all sequences with all *de novo* tandem repeat detection algorithms
defined in the :ref:`configuration file <configure>`:
::

    dTandem_repeats = {}
    for iS in lHIV_Sequence:
        iTandem_repeats = iS.detect(denovo = True)
        dTandem_repeats[iS.id] = iTandem_repeats

    print(dTandem_repeats['sp|Q75006|REV_HV1ET'].repeats[0])


Output the detected tandem repeats.
-----------------------------------

Write a singe repeat_list to .tsv format:
::

    path_to_output_tsv_file = "/my/path/to/the/outputfile.tsv"
    tandem_repeats.write(format = "tsv", file = path_to_output_tsv_file)


Write a singe repeat_list to .pickle format:
::

    path_to_output_pickle_file = "/my/path/to/the/outputfile.pickle"
    tandem_repeats.write(format = "pickle", file = path_to_output_pickle_file)


A repeat_list in pickle format can easily be read in again:
::

    from tandemrepeats.repeat_list import repeat_list
    tandem_repeats = repeat_list.Repeat_list.create(format = "pickle", file = path_to_output_pickle_file)