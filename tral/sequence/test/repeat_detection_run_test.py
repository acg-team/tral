import logging
import logging.config
logging.config.fileConfig("/Users/blieadmin/.tral/logging.ini")

# import pytest

from tral.sequence import repeat_detection_run, sequence

TEST_SEQUENCE_Q9BRR0 = "MARELSESTALDAQSTEDQMELLVIKVEEEEAGFPSSPDLGSEGSRERFRGFRYPEAAGPREALSRLRELCRQWLQPEMHSKEQILELLVLEQFLTILPGNLQSWVREQHPESGEEVVVLLEYLERQLDEPAPQVSGVDQGQELLCCKMALLTPAPGSQSSQFQLMKALLKHESVGSQPLQDRVLQVPVLAHGGCCREDKVVASRLTPESQGLLKVEDVALTLTPEWTQQDSSQGNLCRDEKQENHGSLVSLGDEKQTKSRDLPPAEELPEKEHGKISCHLREDIAQIPTCAEAGEQEGRLQRKQKNATGGRRHICHECGKSFAQSSGLSKHRRIHTGEKPYECEECGKAFIGSSALVIHQRVHTGEKPYECEECGKAFSHSSDLIKHQRTHTGEKPYECDDCGKTFSQSCSLLEHHRIHTGEKPYQCSMCGKAFRRSSHLLRHQRIHTGDKNVQEPEQGEAWKSRMESQLENVETPMSYKCNECERSFTQNTGLIEHQKIHTGEKPYQCNACGKGFTRISYLVQHQRSHVGKNILSQ"

# @pytest.mark.slow


def test_detect_TRUST():

    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    predicted_repeats = repeat_detection_run.run_detector(seq_records=[test_seq], detectors=["TRUST"])[0]['TRUST']
    # Warning: TRUST finds these results ONLY with the BLOSUM50 substitution matrix.
    assert len(predicted_repeats) == 3
    # Warning: TRUST finds these results ONLY with the BLOSUM50 substitution matrix.
    assert predicted_repeats[0].msa == ['HLREDIAQIP---TCAEAGE---QEGRLQR', 'KQKNATGGRR--HICHECGKSFAQSSGLSK', 'HRRIHTGEKP--YECEECGKAFIGSSALVI', 'HQRVHTGEKP--YECEECGKAFSHSSDLIK', 'HQRTHTGEKP--YECDDCGKTFSQSCSLLE', 'HHRIHTGEKP--YQCSMCGKAFRRSSHLLR', 'HQRIHTGDKN--VQEPEQGEAW--KSRM--', 'ESQLENVETPmsYKCNECERSFTQNTGLIE', 'HQKIHTGEKP--YQCNACGKGFTRISYLVQ']


def test_detect_XSTREAM():
    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    predicted_repeats = repeat_detection_run.run_detector(seq_records=[test_seq], detectors=["XSTREAM"])[0]['XSTREAM']
    assert len(predicted_repeats) == 1
    assert predicted_repeats[0].msa == ['ECGKSFAQS-SGLSK-HRRIHTGEKPYECE', 'ECGKAFIGS-SALVI-HQRVHTGEKPYECE', 'ECGKAFSHS-SDL-IKHQRTHTGEKPYECD', 'DCGKTFSQSCSLLEH-H-RIHTGEKPY']


def test_detect_TREKS():

    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    predicted_repeats = repeat_detection_run.run_detector(seq_records=[test_seq], detectors=["T-REKS"])[0]['T-REKS']
    assert len(predicted_repeats) == 1
    assert predicted_repeats[0].msa == ['C---G---KSFAQSSGLSKHRRIHTGEKPYECE-E', 'C---G---KAFIGSSALVIHQRVHTGEKPYECE-E', 'C---G---KAFSHSSDLIKHQRTHTGEKPYECD-D', 'C---G---KTFSQSCSLLEHHRIHTGEKPYQCS-M', 'C---G---KAFRRSSHLLRHQRIHTGDKNVQ-EPE', 'Q---G---EAW--KSRME-SQ-LENVETPMSYK--', 'C---NECERSFTQNTGLIEHQKIHTGEKPYQ----', 'CNACG---KGFTRISYLVQHQRSHVG-KNI-LS--']


def test_detect_HHrepID():

    test_seq = sequence.Sequence(TEST_SEQUENCE_Q9BRR0)
    predicted_repeats = repeat_detection_run.run_detector(seq_records=[test_seq], detectors=["HHrepID"])[0]['HHrepID']
    assert len(predicted_repeats) == 1
    expected_msa = ['------------------IPTCAEAGEQ----',
                    'EGRLQRKQKNATGGRRHICHECGKSFAQ----',
                    'SSGLSKHRRIHTGEKPYECEECGKAFIG----',
                    'SSALVIHQRVHTGEKPYECEECGKAFSH----',
                    'SSDLIKHQRTHTGEKPYECDDCGKTFSQ----',
                    'SCSLLEHHRIHTGEKPYQCSMCGKAFRR----',
                    'SSHLLRHQRIHTGDKNVQEPEQGEAWKSRMES',
                    '------QLENVETPMSYKCNECERSFTQ----',
                    'NTGLIEHQKIHTGEKPYQCNACGKGFTR----',
                    'ISYLVQHQRSHVG-------------------']
    assert predicted_repeats[0].msa == expected_msa
