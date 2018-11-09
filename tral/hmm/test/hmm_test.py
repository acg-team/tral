import os
import pytest

from tral.hmm.hmm import HMM
from tral.repeat import repeat

TEST_SEQUENCE = "GYEDEDEDRPFYALGLGKRPRTYSFGL"

TEST_REPEAT_MSA_DOUBLE = ["AA", "AA"]
TEST_HMM_STATES_DOUBLE = ["N", "M1", "I1", "M2", "I2", "C"]
TEST_HMM_P0_DOUBLE = {i: 0 for i in TEST_HMM_STATES_DOUBLE}


TEST_REPEAT_MSA_SINGLE = ["A", "A", "A"]
TEST_HMM_STATES_SINGLE = ["N", "M1", "I1", "C"]
TEST_HMM_P0_SINGLE = {i: 0 for i in TEST_HMM_STATES_SINGLE}


def test_create_HMM_from_Repeat():

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)

    assert test_hmm.l_effective == 2
    assert set(test_hmm.states) == set(TEST_HMM_STATES_DOUBLE)
    assert test_hmm.p_0 == TEST_HMM_P0_DOUBLE
    #assert test_hmm.p_t == TEST_HMM_P0_DOUBLE

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_SINGLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)

    assert test_hmm.l_effective == 1
    assert test_hmm.states == TEST_HMM_STATES_SINGLE
    assert test_hmm.p_0 == TEST_HMM_P0_SINGLE


def test_create_HMM_from_DNA_Repeat():

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_DOUBLE, sequence_type="DNA")
    #test_parameters = {"hmmbuild":{"hmm_copy_path":"/Users/elkeschaper/Downloads", "hmm_copy_id":"maulwurf"}}
    #test_hmm = HMM.create(input_format = 'repeat', repeat = test_repeat, **test_parameters)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)

    assert set(test_hmm.hmmer['letters']) == set(["A", "C", "G", "T"])
    assert set(test_hmm.alphabet) == set(["A", "C", "G", "T"])

    assert test_hmm.l_effective == 2
    assert set(test_hmm.states) == set(TEST_HMM_STATES_DOUBLE)
    assert test_hmm.p_0 == TEST_HMM_P0_DOUBLE

    for iState, ip_e in test_hmm.p_e.items():
        assert set(list(ip_e.keys())) == set(["A", "C", "G", "T"])
    #assert test_hmm.p_t == TEST_HMM_P0_DOUBLE


def test_hmm_pickle(tmpdir):

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)

    test_pickle = tmpdir.join("test.pickle")
    test_hmm.write(test_pickle, 'pickle')
    test_hmm_new = HMM.create(input_format='pickle', file=test_pickle)

    assert test_hmm.hmmer == test_hmm_new.hmmer
    assert test_hmm.alphabet == test_hmm_new.alphabet

    if os.path.exists(test_pickle):
        os.remove(test_pickle)
