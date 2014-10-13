import os
import pytest

from tandemrepeats.hmm.hmm import HMM
from tandemrepeats.repeat import repeat

TEST_SEQUENCE = "GYEDEDEDRPFYALGLGKRPRTYSFGL"

TEST_REPEAT_MSA_DOUBLE = ["AA","AA"]
TEST_HMM_STATES_DOUBLE = ["N", "M1", "I1","M2","I2","C"]
TEST_HMM_P0_DOUBLE = {i:0 for i in TEST_HMM_STATES_DOUBLE}


TEST_REPEAT_MSA_SINGLE = ["A","A","A"]
TEST_HMM_STATES_SINGLE = ["N", "M1", "I1","C"]
TEST_HMM_P0_SINGLE = {i:0 for i in TEST_HMM_STATES_SINGLE}

# Test file names
TEST_FILE_WITH_ID = 'carcinustatin.hmm'

@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'hmm', 'test')


def test_create_HMM_from_Repeat():


    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(format = 'repeat', repeat = test_repeat)

    assert test_hmm.lD == 2
    assert set(test_hmm.states) == set(TEST_HMM_STATES_DOUBLE)
    assert test_hmm.p_0 == TEST_HMM_P0_DOUBLE
    #assert test_hmm.p_t == TEST_HMM_P0_DOUBLE


    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA_SINGLE)
    test_hmm = HMM.create(format = 'repeat', repeat = test_repeat)

    assert test_hmm.lD == 1
    assert test_hmm.states == TEST_HMM_STATES_SINGLE
    assert test_hmm.p_0 == TEST_HMM_P0_SINGLE


def test_hmm_pickle():

    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(format = 'repeat', repeat = test_repeat)

    test_pickle = os.path.join(path(), "test.pickle")
    test_hmm.write(test_pickle, 'pickle')
    test_hmm_new = HMM.create(format = 'pickle', file = test_pickle)

    assert test_hmm.hmmer == test_hmm_new.hmmer
    assert test_hmm.alphabet == test_hmm_new.alphabet

    if os.path.exists(test_pickle):
        os.remove(test_pickle)