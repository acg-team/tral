import collections
import logging
import os
import pytest

from tandemrepeats.repeat_list import repeat_list as rl
from tandemrepeats.repeat import repeat

TEST_REPEATS = [["AA","AA"],["AAA","AAA"],["AAAA","AAAA"], ["AAA-","AAAA"]]
TEST_SCORE = "phylo_gap01"
TEST_SCORE_VALUE_LIST = [0.0, 0.5, 1.0,  1.0]
TEST_BEGIN_LIST = [6,10,10,10]

def test_create_Repeat_list_from_Repeats():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    test_repeat_list = rl.Repeat_list(repeats = test_repeats)

    assert len(test_repeat_list.repeats) == 4
    for i,j in zip (TEST_REPEATS, test_repeat_list.repeats):
        assert i == j.msa


def test_pairwise_overlap():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j

    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[:2]) == False
    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[1:3]) == False
    assert rl.two_repeats_overlap("common_ancestry", *test_repeats[2:]) == False
    assert rl.two_repeats_overlap("shared_char", *test_repeats[:2]) == False
    assert rl.two_repeats_overlap("shared_char", *test_repeats[1:3]) == True


def test_cluster():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j

    test_repeat_list = rl.Repeat_list(repeats = test_repeats)
    test_repeat_list.cluster("common_ancestry")

    # Check whether both lists include exactly the same elements.
    for i in [{0}, {1,3}, {2}]:
        assert i in test_repeat_list.dCluster["common_ancestry"]
    assert len(test_repeat_list.dCluster["common_ancestry"]) == 3

    test_repeat_list.cluster("shared_char")

    # Check whether both lists include exactly the same elements.
    for i in [{0}, {1,2,3}]:
        assert i in test_repeat_list.dCluster["shared_char"]
    assert len(test_repeat_list.dCluster["shared_char"]) == 2


def test_filter_pValue():

    #test_repeats = [repeat.Repeat(msa = i, scoreslist = ["phylo_gap01"], calc_score = True, calc_pValue = True) for i in TEST_REPEATS]
    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_SCORE_VALUE_LIST):
        i.dPValue = {}
        i.dPValue[TEST_SCORE] = j

    test_repeat_list = rl.Repeat_list(repeats = test_repeats)

    test_repeat_list_filtered = test_repeat_list.filter("pValue", TEST_SCORE, 0.1)
    assert len(test_repeat_list_filtered.repeats) == 1


def test_filter_cluster_based():

    test_repeats = [repeat.Repeat(msa = i) for i in TEST_REPEATS]
    for i,j in zip(test_repeats, TEST_SCORE_VALUE_LIST):
        i.dPValue = {}
        i.dPValue[TEST_SCORE] = j
    for i,j in zip(test_repeats, TEST_BEGIN_LIST):
        i.begin = j


    test_repeat_list = rl.Repeat_list(repeats = test_repeats)
    test_repeat_list.filter("pValue", TEST_SCORE, 0.1)
    test_repeat_list_filtered = test_repeat_list.filter("none_overlapping", ["common_ancestry"], {"pValue":TEST_SCORE, "divergence":TEST_SCORE})
    assert len(test_repeat_list_filtered.repeats) == 3
    for i in test_repeats[:3]:
        assert i in test_repeat_list_filtered.repeats