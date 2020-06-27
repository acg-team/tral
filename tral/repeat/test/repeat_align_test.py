# import pytest

from tral.repeat import repeat_align

TEST_MSA = ["MGKGYL---------------------------------------ALCSYNCKEA-INILSHLPSHHYN",
            "TG--------------------------------------------------------------WVLCQ",
            "IGRAYF---------------------------------------ELSEYMQAER-IFSEVRRIENYRV",
            "EGMEIYSTTLWHLQK------------------------------DVALSVLSKDLTDMDKNSPEAWCA",
            "AGNCFS---------LQREH-------------------------DIAIKFFQRA-IQVDPNYAYAYTL",
            "LGHEFV--------------LTEEL--------------------DKALACFRNA-IRVNPRHYNAWYG",
            "LGMIYY-------------------KQEKF---------------SLAEMHFQKA-LDINPQSSVLLCH",
            "IGVVQH------------------------ALKKS----------EKALDTLNKA-IVIDPKNPLCKFH",
            "RASVLF-----------------------------ANEKY-----KSALQELEEL-KQIVPKESLVYFL",
            "IGKVYK----------------------------------KLGQTHLALMNFSWA-MDLDPKGAN----"]
TEST_REALIGNED_MSA = ['MGKGY-----------LALCSYNCKEAI-NILSHLPSHHYN',
                      'TGWVL----------------------------------CQ',
                      'IGRAYFELS----EYMQAERIFS--EVR-RIEN-----YRV',
                      'EGMEIYSTTLWHLQKDVALSVLS--KDLTDMDKNSPEAWCA',
                      'AGNCFSLQR----EHDIAIKFFQ--RAI-QVDPNYAYAYTL',
                      'LGHEFVLTE----ELDKALACFR--NAI-RVNPRHYNAWYG',
                      'LGMIYYKQE----KFSLAEMHFQ--KAL-DINPQSSVLLCH',
                      'IGVVQHALK----KSEKALDTLN--KAI-VIDPKNPLCKFH',
                      'RASVLFANE----KYKSALQELE--ELK-QIVPKESLVYFL',
                      'IGKVYKKLG----QTHLALMNFS--WAM-DLDPKGA----N']


def test_repeat_alignment():
    ''' Test a realignment using Mafft's ginsi.'''

    REALIGNED_MSA = repeat_align.realign_repeat(my_msa=TEST_MSA)

    assert isinstance(REALIGNED_MSA, list)

    # Attention: If ginsi changes, the alignment will also change...
    assert REALIGNED_MSA == TEST_REALIGNED_MSA
