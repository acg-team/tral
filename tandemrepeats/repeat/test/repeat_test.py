import numpy as np
import pytest

from tandemrepeats.repeat import repeat
from tandemrepeats.repeat import repeat_io


TEST_MSA_O = ['OCC', 'OOO']
TEST_MSA_K = ['KCC', 'KKK']
TEST_SCORE = "phylo"

# defaultdict(<class 'int'>, {'pSim': 0.66666666669999997, 'parsimony': 0.66666666669999997, 'entropy': 0.66666666666666663, 'phylo': 0.11368675605567802})
# pValue 'phylo': 0.3821

notfixed = pytest.mark.notfixed

def test_standardize_amino_acids():

    assert repeat.standardize("ABDEF-G") == "ADDEF-G"


def test_repeat_ambiguous():

    myTR_O = repeat.Repeat(msa = TEST_MSA_O)
    myTR_K = repeat.Repeat(msa = TEST_MSA_K)

    assert myTR_O.msaTD_standard_aa == myTR_K.msaTD
    assert myTR_O.msaTD_standard_aa == myTR_K.msaTD_standard_aa

    assert myTR_O.score(TEST_SCORE) == myTR_K.score(TEST_SCORE)
    assert myTR_O.divergence(TEST_SCORE) == myTR_K.divergence(TEST_SCORE)
    assert myTR_O.pValue(TEST_SCORE) == myTR_K.pValue(TEST_SCORE)

    assert myTR_O.divergence(TEST_SCORE) == 2.5986328125
    assert myTR_K.pValue(TEST_SCORE) == 0.3821