import os
import pytest

from tral.repeat import repeat


TEST_MSA_O = ['OCC', 'OOO']
TEST_MSA_K = ['KCC', 'KKK']
TEST_SCORE = "phylo_gap01"

# defaultdict(<class 'int'>, {'pSim': 0.66666666669999997, 'parsimony': 0.66666666669999997,
#                             'entropy': 0.66666666666666663, 'phylo': 0.11368675605567802})
# pvalue 'phylo': 0.3821


@pytest.mark.no_external_software_required
def test_standardize_amino_acids():

    assert repeat.standardize("ABDEF-G", "AA") == "ADDEF-G"


@pytest.mark.no_external_software_required
def test_repeat_ambiguous():

    myTR_O = repeat.Repeat(msa=TEST_MSA_O)
    myTR_K = repeat.Repeat(msa=TEST_MSA_K)

    assert myTR_O.msaTD_standard_aa == myTR_K.msaTD
    assert myTR_O.msaTD_standard_aa == myTR_K.msaTD_standard_aa

    assert myTR_O.score(TEST_SCORE) == myTR_K.score(TEST_SCORE)
    assert myTR_O.divergence(TEST_SCORE) == myTR_K.divergence(TEST_SCORE)
    assert myTR_O.pvalue(TEST_SCORE) == myTR_K.pvalue(TEST_SCORE)

    assert myTR_O.divergence(TEST_SCORE) == 2.095947265625
    assert myTR_K.pvalue(TEST_SCORE) == 0.3507


@pytest.mark.no_external_software_required
def test_repeat_pickle(tmpdir):

    myTR_O = repeat.Repeat(msa=TEST_MSA_O)

    test_pickle = os.path.join(tmpdir,"test.pickle")
    myTR_O.write(test_pickle, 'pickle')
    myTR_O_new = repeat.Repeat.create(test_pickle, 'pickle')

    assert myTR_O.msa == myTR_O_new.msa
    assert myTR_O.sequence_type == myTR_O_new.sequence_type
    assert myTR_O.text == myTR_O_new.text

    if os.path.exists(test_pickle):
        os.remove(test_pickle)

