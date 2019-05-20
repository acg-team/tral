import os
import pytest

from tral.hmm.hmm import HMM
from tral.sequence import sequence
from tral.repeat import repeat
from tral.repeat_list import repeat_list


TEST_SEQUENCE = "FFAAAAAAFF"
TEST_SEQUENCE_CARC = "FFAGPYAYGLAGPYAYGLFF"
# Zinc finger protein Q9BRR0
TEST_SEQUENCE_Q9BRR0 = "MARELSESTALDAQSTEDQMELLVIKVEEEEAGFPSSPDLGSEGSRERFRGFRYPEAAGPREALSRLRELCRQWLQPEMHSKEQILELLVLEQFLTILPGNLQSWVREQHPESGEEVVVLLEYLERQLDEPAPQVSGVDQGQELLCCKMALLTPAPGSQSSQFQLMKALLKHESVGSQPLQDRVLQVPVLAHGGCCREDKVVASRLTPESQGLLKVEDVALTLTPEWTQQDSSQGNLCRDEKQENHGSLVSLGDEKQTKSRDLPPAEELPEKEHGKISCHLREDIAQIPTCAEAGEQEGRLQRKQKNATGGRRHICHECGKSFAQSSGLSKHRRIHTGEKPYECEECGKAFIGSSALVIHQRVHTGEKPYECEECGKAFSHSSDLIKHQRTHTGEKPYECDDCGKTFSQSCSLLEHHRIHTGEKPYQCSMCGKAFRRSSHLLRHQRIHTGDKNVQEPEQGEAWKSRMESQLENVETPMSYKCNECERSFTQNTGLIEHQKIHTGEKPYQCNACGKGFTRISYLVQHQRSHVGKNILSQ"


TEST_REPEAT_MSA_SINGLE = ["A", "A", "A"]
TEST_RESULT_REPEAT_MSA_SINGLE = ["A", "A", "A", "A", "A", "A"]

TEST_REPEAT_MSA_DOUBLE = ["AA", "AA"]
TEST_RESULT_REPEAT_MSA_DOUBLE = ["AA", "AA", "AA"]

TEST_RESULT_REPEAT_MSA_LONG = [10 * "A", 10 * "A"]
TEST_SEQUENCE_A = 18 * "A"

TEST_RESULT_REPEAT_MSA_SUPER_LONG = [100 * "A", 100 * "A"]
TEST_SEQUENCE_SUPER_LONG_A = 200 * "A"


# Test file names
TEST_FILE_WITH_ID = 'carcinustatin.hmm'

TEST_SEQUENCE_TAG = 'test_sequence_tag'


@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.mark.no_external_software_required
def test_initialise_sequence():
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    assert test_seq.seq == TEST_SEQUENCE


@pytest.mark.no_external_software_required
def test_detect_repeats_with_hmm(path):
    test_hmm = HMM.create(input_format='hmmer', file=os.path.join(path, "..", "..", "hmm", "test", TEST_FILE_WITH_ID))
    test_seq = sequence.Sequence(TEST_SEQUENCE_CARC)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 1


def test_detect_repeats_with_repeat():

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 1
    assert test_optimized_repeat.repeats[0].msa == TEST_RESULT_REPEAT_MSA_DOUBLE

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_SINGLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 1
    assert test_optimized_repeat.repeats[0].msa == TEST_RESULT_REPEAT_MSA_SINGLE


@pytest.mark.no_external_software_required
def test_too_big_hmms():

    test_repeat = repeat.Repeat(msa=TEST_RESULT_REPEAT_MSA_LONG)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)
    test_seq = sequence.Sequence(TEST_SEQUENCE_A)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 0

    test_repeat = repeat.Repeat(msa=TEST_RESULT_REPEAT_MSA_SUPER_LONG)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)
    test_seq = sequence.Sequence(TEST_SEQUENCE_SUPER_LONG_A)
    test_optimized_repeat = test_seq.detect([test_hmm])
    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 0


def test_detect_repeats_denovo():

    test_parameters = {"detection": {"detectors": ["TRUST"]}}

    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    test_optimized_repeat = test_seq.detect(denovo=True, **test_parameters)

    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert len(test_optimized_repeat.repeats) == 3


@pytest.mark.no_external_software_required
def test_sequence_pickle(tmpdir):

    test_seq = sequence.Sequence(TEST_SEQUENCE)

    test_pickle = tmpdir.join("test.pickle")
    test_seq.write(test_pickle, 'pickle')
    test_seq_new = sequence.Sequence.create(test_pickle, 'pickle')

    assert test_seq.seq == test_seq_new.seq

    test_repeat = repeat.Repeat(msa=TEST_REPEAT_MSA_DOUBLE)
    test_hmm = HMM.create(input_format='repeat', repeat=test_repeat)
    test_optimized_repeat = test_seq.detect([test_hmm])
    test_seq.set_repeatlist(test_optimized_repeat, TEST_SEQUENCE_TAG)

    assert type(test_optimized_repeat) == repeat_list.RepeatList
    assert list(test_seq.d_repeatlist.keys()) == [TEST_SEQUENCE_TAG]
    assert type(test_seq.d_repeatlist[TEST_SEQUENCE_TAG]) == repeat_list.RepeatList
    assert test_seq.d_repeatlist[TEST_SEQUENCE_TAG].repeats

    test_retrieved_repeatlist = test_seq.get_repeatlist(TEST_SEQUENCE_TAG)
    assert test_retrieved_repeatlist == test_optimized_repeat

    test_seq.write(test_pickle, 'pickle')
    test_seq_new = sequence.Sequence.create(test_pickle, 'pickle')

    assert test_seq.d_repeatlist.keys() == test_seq_new.d_repeatlist.keys()
    assert test_seq.d_repeatlist[TEST_SEQUENCE_TAG].repeats[0].msa == test_seq_new.d_repeatlist[TEST_SEQUENCE_TAG].repeats[0].msa

    if os.path.exists(test_pickle):
        os.remove(test_pickle)
