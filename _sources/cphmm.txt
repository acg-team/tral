.. _cphmm:

Annotate tandem repeats from sequence domain models.
====================================================

Here you learn how to annotate sequences with your favourite sequence profile model.
The sequence profile model is transformed into a circular profile HMM, as described in
a :ref:`recent MBE publication <publications>`. Next, the sequence
is searched for tandem repeats using the circular profile HMM.



Read in your sequence profile model.
------------------------------------

::

    import os
    from tandemrepeats.hmm import hmm
    from tandemrepeats.paths import *

    fPfam_profile_hmm = os.path.join(PACKAGE_DIRECTORY,"test","zf-CCHC.hmm")

    circular_profile_HMM_zfCCHC = hmm.HMM.create(format = 'hmmer', file = fPfam_profile_hmm)

    print(circular_profile_HMM_zfCCHC)



Read in your sequences.
-----------------------

::

    from tandemrepeats.sequence import sequence

    human_HCFC1_fasta = os.path.join(path_to_tandemrepeats,"test","P62633.fasta")
    human_HCFC1_sequence = sequence.Sequence.create(file = human_HCFC1_fasta, format = 'fasta')[0]




Annotate tandem repeats with the circular profile HMM.
------------------------------------------------------

::

    tandem_repeats = human_HCFC1_sequence.detect(lHMM = [circular_profile_HMM_zfCCHC])

    print(tandem_repeats.repeats[0])



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