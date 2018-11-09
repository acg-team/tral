import collections
import os
import pytest

from tral.repeat_list import repeat_list as rl
from tral.repeat import repeat
from tral.repeat_list import repeat_list_io as rl_io
from tral.sequence import sequence

TEST_REPEATS = [["AA", "AA"], ["AAA", "AAA"]]
TEST_SEQUENCE = "MAAAAKAAAAAAL"
TEST_SCORE = "phylo_gap01"
TEST_SCORE_VALUE_LIST = [0.0, 0.5]
TEST_BEGIN_LIST = [6, 10]
# The resulting string should contain the following data, however perhaps in a different order:
TEST_TSV = "msa_original\tbegin\tn_effective\tl_effective\tsequence_length\tpvalue\nAA,AA\t2\t2.0\t2\t4\tNone\nAAA,AAA\t7\t2.0\t3\t6\tNone"


@pytest.mark.no_external_software_required
def test_serialize_repeat_list_tsv():

    test_repeats = [repeat.Repeat(msa=i) for i in TEST_REPEATS]
    test_seq = sequence.Sequence(TEST_SEQUENCE)
    for i in test_repeats:
        test_seq.repeat_in_sequence(i)
    test_repeat_list = rl.RepeatList(repeats=test_repeats)

    tsv = rl_io.serialize_repeat_list_tsv(test_repeat_list)

    assert type(tsv) == str
