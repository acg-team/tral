import logging
import os
import pytest

from tandemrepeats.hmm.hmm import HMM
from tandemrepeats.sequence import sequence
from tandemrepeats.repeat import repeat
from tandemrepeats.repeat_list import repeat_list


TEST_SEQUENCE = "FFAAAAAAFF"
# Zinc finger protein Q9BRR0
TEST_SEQUENCE_Q9BRR0 = "MARELSESTALDAQSTEDQMELLVIKVEEEEAGFPSSPDLGSEGSRERFRGFRYPEAAGPREALSRLRELCRQWLQPEMHSKEQILELLVLEQFLTILPGNLQSWVREQHPESGEEVVVLLEYLERQLDEPAPQVSGVDQGQELLCCKMALLTPAPGSQSSQFQLMKALLKHESVGSQPLQDRVLQVPVLAHGGCCREDKVVASRLTPESQGLLKVEDVALTLTPEWTQQDSSQGNLCRDEKQENHGSLVSLGDEKQTKSRDLPPAEELPEKEHGKISCHLREDIAQIPTCAEAGEQEGRLQRKQKNATGGRRHICHECGKSFAQSSGLSKHRRIHTGEKPYECEECGKAFIGSSALVIHQRVHTGEKPYECEECGKAFSHSSDLIKHQRTHTGEKPYECDDCGKTFSQSCSLLEHHRIHTGEKPYQCSMCGKAFRRSSHLLRHQRIHTGDKNVQEPEQGEAWKSRMESQLENVETPMSYKCNECERSFTQNTGLIEHQKIHTGEKPYQCNACGKGFTRISYLVQHQRSHVGKNILSQ"


TEST_REPEAT_MSA_SINGLE = ["A","A","A"]
TEST_RESULT_REPEAT_MSA_SINGLE = ["A","A","A","A","A","A"]

TEST_REPEAT_MSA_DOUBLE = ["AA","AA"]
TEST_RESULT_REPEAT_MSA_DOUBLE = ["AA","AA","AA"]

#Test file names
TEST_FILE_WITH_ID = 'carcinustatin.hmm'

notfixed = pytest.mark.notfixed

@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'hmm', 'test')

def test_initialise_sequence():
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    assert test_seq.seq == TEST_SEQUENCE

def test_detect_repeats_with_hmm():
    test_hmm = HMM.create(hmmer_file = os.path.join(path(), TEST_FILE_WITH_ID))
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    test_optimized_repeat = test_seq.detect([test_hmm])

#@notfixed
def test_detect_repeats_with_repeat():

    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(repeat = test_repeat)
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.Repeat_list
    assert len(test_optimized_repeat.repeats) == 1
    assert test_optimized_repeat.repeats[0].msa == TEST_RESULT_REPEAT_MSA_DOUBLE

    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA_SINGLE)
    test_hmm = HMM.create(repeat = test_repeat)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.Repeat_list
    assert len(test_optimized_repeat.repeats) == 1
    assert test_optimized_repeat.repeats[0].msa == TEST_RESULT_REPEAT_MSA_SINGLE

#@notfixed
def test_detect_repeats_denovo():

    test_parameters = {"detection": {"lFinders": ["TRUST"]}}

    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    test_optimized_repeat = test_seq.detect(denovo = True, **test_parameters)

    assert type(test_optimized_repeat) == repeat_list.Repeat_list
    assert len(test_optimized_repeat.repeats) == 3