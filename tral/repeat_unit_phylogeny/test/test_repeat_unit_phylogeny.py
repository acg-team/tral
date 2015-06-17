import os
import pytest

from tral.repeat_unit_phylogeny import repeat_unit_phylogeny
from tral.repeat_unit_phylogeny.repeat_unit_phylogeny import RepeatUnitPhylogeny
from tral.repeat import repeat
from tral.hmm import hmm, hmm_viterbi

# In case all repeat units start from the same index
TEST_REPEAT_MSA_ORTHOLOG1 = ["AA","AA"]
TEST_REPEAT_MSA_ORTHOLOG2 = ["A","A","A"]
TEST_REPEAT_MSA_ORTHOLOG3 = ["A","A","A"]

TEST_REPEAT_MSA = ["AC","AC"]
TEST_SEQUENCE1 = "GYEDEACACACLPRY"
TEST_SEQUENCE2 = "DEDRATACACLGLG"
TEST_SEQUENCE3 = "FYCACACACACACALGLGKRP"
TEST_IDS = ["S1", "S2", "S3"]

# In case repeat units were detected with the Viterbi algorithm and do not
# from the same index
TEST_SEQUENCES = ""
TEST_VITERBI_PATHS1 =
TEST_VITERBI_PATHS2 =
TEST_VITERBI_PATHS3 =


@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'repeat_unit_phylogeny', 'test')


@pytest.fixture
def orthologs():
    """Return repeats
    """
    return [repeat.Repeat(msa = i) for i in [TEST_REPEAT_MSA_ORTHOLOG1, TEST_REPEAT_MSA_ORTHOLOG2, TEST_REPEAT_MSA_ORTHOLOG3]]

@notfixed
def test_create_RepeatUnitPhylogeny_from_single_Repeat():

    test_repeat = orthologs()[0]
    test_phylo = RepeatUnitPhylogeny.create(input_type = 'repeat', repeats = test_repeat)
    assert 1 == 1

@notfixed
def test_create_RepeatUnitPhylogeny_from_double_Repeat():

    test_repeats = orthologs()[0:1]
    test_phylo = RepeatUnitPhylogeny.create(input_type = 'repeat', repeats = test_repeats)
    assert 1 == 1

@notfixed
def test_check_Repeat_characteristics():

    assert 1 == 1


def pair_wise_phylogenies():

    # Retrieve Viterbi paths
    test_repeat = repeat.Repeat(msa = TEST_REPEAT_MSA)
    test_hmm = hmm.HMM.create(input_format = 'repeat', repeat = test_repeat)
    test_viterbi_path1 = hmm_viterbi.viterbi(test_hmm, TEST_SEQUENCE1)
    test_viterbi_path2 = hmm_viterbi.viterbi(test_hmm, TEST_SEQUENCE2)
    test_viterbi_path3 = hmm_viterbi.viterbi(test_hmm, TEST_SEQUENCE3)

    viterbi_paths = [test_viterbi_path1, test_viterbi_path2, test_viterbi_path3]
    test_sequences = [TEST_SEQUENCE1, TEST_SEQUENCE2, TEST_SEQUENCE3]
    lD = 2

    test_d_sequences = {iID: iS for iID, iS in zip(TEST_IDS, test_sequences)}
    test_d_viterbi_paths = {iID: iV for iID, iV in zip(TEST_IDS, viterbi_paths)}

    test_d_phylo = repeat_unit_phylogeny.pairwise_optimal_repeat_unit_phylogenies(test_d_sequences, test_d_viterbi_paths, lD)

    assert ('S2', 'S3') in test_d_phylo
    assert isinstance(test_d_phylo[('S2', 'S3')].phylogeny, ete3.coretype.tree.TreeNode)
