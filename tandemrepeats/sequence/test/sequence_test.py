import logging
import os
import pytest

from tandemrepeats.hmm.hmm import HMM
from tandemrepeats.sequence import sequence

TEST_SEQUENCE = "GYEDEDEDRPFYALGLGKRPRTYSFGL"


# Test file names
TEST_FILE_WITH_ID = 'carcinustatin.hmm'

@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'hmm', 'test')


def test_initialise_sequence(sequence):

    test_seq = seq.Seq(sequence)


def test_detect_repeats_with_hmm(path):
    test_hmm = HMM.create(hmmer_file = os.path.join(path, TEST_FILE_WITH_ID))
    test_seq = seq.Seq(sequence)
    test_optimized_repeat = test_seq.detect(test_hmm)

    #assert
