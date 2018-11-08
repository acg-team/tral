import collections
import os
import pytest

from tral.repeat_list import repeat_list as rl
from tral.repeat import repeat
from tral.sequence import sequence
from tral.paths import PACKAGE_DIRECTORY

TEST_REPEATS = [["AA","AA"],["AAA","AAA"],["AAAA","AAAA"], ["AAA-","AAAA"]]
TEST_SCORE = "phylo_gap01"
TEST_SCORE_VALUE_LIST = [0.0, 0.5, 1.0,  1.0]
TEST_BEGIN_LIST = [6,10,10,10]
TEST_SEQUENCE = "MAAAAKAAAAAAL"

# The resulting string should contain the following data, however perhaps in a different order:
TEST_TSV = "msa_original\tbegin\tn_effective\tl_effective\tsequence_length\tpvalue\nAA,AA\t2\t2.0\t2\t4\tNone\nAAA,AAA\t7\t2.0\t3\t6\tNone"

@pytest.mark.no_external_software_required
def test_create_repeat_list_from_repeats():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    test_repeat_list = rl.RepeatList(repeats = test_repeats)

    assert len(test_repeat_list.repeats) == 4
    for i,j in zip (TEST_REPEATS, test_repeat_list.repeats):
        assert i == j.msa


@pytest.mark.no_external_software_required
def test_repeat_list_pickle(tmpdir):
    """Tests saving repeats as a pickle.
    
    Args:
        - tmpdir (py.path.local) Fixture to generate a temporary directory
    """
    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    test_repeat_list = rl.RepeatList(repeats = test_repeats)

    #test_pickle = os.path.join(path, "test.pickle")
    test_pickle = tmpdir.join("test.pickle")
    test_repeat_list.write('pickle', test_pickle)
    test_repeat_list_new = repeat.Repeat.create(test_pickle, 'pickle')

    assert len(test_repeat_list.repeats) == len(test_repeat_list_new.repeats)
    assert test_repeat_list.repeats[0].msa == test_repeat_list_new.repeats[0].msa

    if os.path.exists(test_pickle):
        os.remove(test_pickle)


def test_serialize_repeat_list_tsv():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS[:2]]
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    for i in test_repeats:
        test_seq.repeat_in_sequence(i)
    test_repeat_list = rl.RepeatList(repeats = test_repeats)

    tsv = test_repeat_list.write("tsv", return_string = True)

    assert type(tsv) == str


@pytest.mark.no_external_software_required
def test_pairwise_overlap():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j

    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[:2]) == False
    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[1:3]) == False
    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[2:]) == False
    assert rl.two_repeats_overlap("shared_char", *test_repeats[:2]) == False
    assert rl.two_repeats_overlap("shared_char", *test_repeats[1:3]) == True


@pytest.mark.no_external_software_required
def test_cluster():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j

    test_repeat_list = rl.RepeatList(repeats = test_repeats)
    test_repeat_list.cluster("common_ancestry")

    # Check whether both lists include exactly the same elements.
    for i in [{0}, {1,3}, {2}]:
        assert i in test_repeat_list.d_cluster["common_ancestry"]
    assert len(test_repeat_list.d_cluster["common_ancestry"]) == 3

    test_repeat_list.cluster("shared_char")

    # Check whether both lists include exactly the same elements.
    for i in [{0}, {1,2,3}]:
        assert i in test_repeat_list.d_cluster["shared_char"]
    assert len(test_repeat_list.d_cluster["shared_char"]) == 2


@pytest.mark.no_external_software_required
def test_filter_pvalue():

    #test_repeats = [repeat.Repeat(msa = i, scoreslist = ["phylo_gap01"], calc_score = True, calc_pvalue = True) for i in TEST_REPEATS]
    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_SCORE_VALUE_LIST):
        i.d_pvalue = {}
        i.d_pvalue[TEST_SCORE] = j

    test_repeat_list = rl.RepeatList(repeats = test_repeats)

    test_repeat_list_filtered = test_repeat_list.filter("pvalue", TEST_SCORE, 0.1)
    assert len(test_repeat_list_filtered.repeats) == 1


@pytest.mark.no_external_software_required
def test_filter_cluster_based():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_SCORE_VALUE_LIST):
        i.d_pvalue = {}
        i.d_pvalue[TEST_SCORE] = j
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j


    test_repeat_list = rl.RepeatList(repeats = test_repeats)
    test_repeat_list.filter("pvalue", TEST_SCORE, 0.1)
    test_repeat_list_filtered = test_repeat_list.filter("none_overlapping", ("common_ancestry", None), [("pvalue", TEST_SCORE), ("divergence", TEST_SCORE)])
    assert len(test_repeat_list_filtered.repeats) == 3
    for i in test_repeats[:3]:
        assert i in test_repeat_list_filtered.repeats