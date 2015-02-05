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
    from tandemrepeats.paths import *

    fHIV_proteome = os.path.join(PACKAGE_DIRECTORY,"test","HIV-1_388796.faa")
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

    for iSequence in lHIV_Sequence:
        iTandem_repeats = iSequence.detect(denovo = True)
        iSequence.set_repeat_list(iTandem_repeats, "denovo")

    print(lHIV_Sequence[0].dRepeat_list['denovo'].repeats[0])


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


Save multiple sequence together with tandem repeat annotations as:
::

    import pickle
    path_to_output_pickle_file = "/my/path/to/the/outputfile.pickle"
    with open(path_to_output_pickle_file, 'wb') as fh:
        pickle.dump(lHIV_Sequence, fh)
