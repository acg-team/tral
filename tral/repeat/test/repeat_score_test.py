import pytest

from tral.repeat import repeat
from tral.repeat import repeat_io


TEST_MSA = ['HPFGFV-------------AVPTKNP-DGTMNLMNWECAIPGKKGTPWEGGLFKLRMLFKDDYPS---SPPKCKFEPPLFHPNV', 'YPSGTVCLsileedkdwrpAITIKQIlLGIQELLN-E---PNIQ-DPAQAEAYTIYCQNRVEYEKrvrAQAK-KFAP-------']
TEST_BEGIN = 0


@pytest.mark.notfixed
@pytest.mark.no_external_software_required
def test_repeat_score():

    # MAKE THIS TEST A TEST!
    Q, eqFreq, alphabet = repeat_io.loadModel()
    indelRatePerSite = 0.001

    myTR = repeat.Repeat(msa=TEST_MSA, begin=TEST_BEGIN)
    myTR.deleteInsertionColumns()

    print(loglikelihood_gaps_starphylogeny_zipfian(t=1, tandem_repeat=myTR))
    print(optimisation(function=loglikelihood_substitution, args=[Q, eqFreq, alphabet, myTR]))
    print(optimisation(function=loglikelihood_substitutions_gaps, args=[[Q, eqFreq, alphabet, myTR], [myTR, indelRatePerSite]]))

    assert 1 == 2
