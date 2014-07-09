import logging
import os
import pytest

from tandemrepeats.sequence import repeat_detection_run, sequence
from tandemrepeats.hmm import hmm

TEST_FAA_FILE = "HIV-1_388796.faa"
TEST_HMM_FILE = "zf-CCHC.hmm"

TEST_DENOVO_PARAMETERS = {"detection": {"lFinders": ["XSTREAM", "TREKS", "TRUST"]}}
TEST_SCORE = "phylo_gap01"

TEST_RESULT_SEQ1 = [3, [["LFNSTKLE","LFNSST-N"], ["GDII", "GDIR"], ["FLG","FLG"]], 3, [["FNCGG-EF","FYCNTSNL","FNSTKLEL","FNSST-NL"], ["GDII", "GDIR"], ["FLG","FLG"]]]

notfixed = pytest.mark.notfixed

@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'test')

#@notfixed
def test_pipeline():

    test_lSeq = sequence.Sequence.read(os.path.join(path(), TEST_FAA_FILE))

    for i,iSeq in enumerate(test_lSeq[:2]):
        test_repeat_denovo = iSeq.detect(denovo = True, **TEST_DENOVO_PARAMETERS)
        test_repeat_denovo_filtered = test_repeat_denovo.filter("pValue", TEST_SCORE, 0.05)
        test_repeat_denovo_hmm = [hmm.HMM.create(repeat = iTR) for iTR in test_repeat_denovo_filtered.repeats]
        test_repeat_denovo_remastered = iSeq.detect(lHMM = test_repeat_denovo_hmm)
        if i == 0:
            assert TEST_RESULT_SEQ1[0] == len(test_repeat_denovo.repeats)
            for k,l in zip(TEST_RESULT_SEQ1[1], test_repeat_denovo.repeats):
                assert k == l.msa
            assert TEST_RESULT_SEQ1[2] == len(test_repeat_denovo_filtered.repeats)
            for k,l in zip(TEST_RESULT_SEQ1[3], test_repeat_denovo_remastered.repeats):
                assert k == l.msa
            # CURRENTLY, ONLY THE FIRST SEQ IS TESTED.
            break

    # CONTINUE WITH SEQUENCE 2
    test_repeat_denovo_remastered_filtered = test_repeat_denovo_remastered.filter("pValue", TEST_SCORE, 0.05)

    test_model_PFAM = hmm.HMM("../../../hmm.faa")
    test_repeat_PFAM = iSeq.detect(lHMM = [test_model_PFAM])
    repeat_PFAM_filtered = repeat_PFAM.filter("pValue", TEST_SCORE, 0.05)

    # Filter test_repeat_denovo_remastered_filtered for none overlapping

#     repeat_optimised_filtered = repeat_optimised.filter(opts)
#
#     final_repeat = repeat_optimised_filtered + repeat_PFAM_filtered
#     final_repeat.write(filehandle, opts)
#
#     model_from_repeat = HMM.create(repeat, opts)
#     repeat_optimised = seq.detect(model_from_repeat, opts)