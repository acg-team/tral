from tral.hmm import hmm
from tral.sequence import sequence
import logging
import logging.config
import os
import pytest

from tral.paths import config_file

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')


TEST_FAA_FILE_MBE_2014 = "P51610.fasta"
TEST_HMM_FILES_MBE_2014 = ["Kelch_1.hmm", "Kelch_2.hmm"]
TEST_SCORE_MBE_2014 = "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions"

# TEST RESULT MBE 2014
# all[all$Ensembl_Protein_ID == 'ENSP00000309555',]
# TR_ID lHMM lRepeat_Region lRepeat_Region_Mean lRepeat_Region_Max lMean  n        n_effective detection_ID
# 19571_1   47            306           257.00353                364  50.0  5  4.061224      PF01344
# 19571_0   49            313           274.15901                354  48.8  6  5.265306      PF07646
# 19571_3   28            426            69.73498                530  29.0 14 13.275862      xstream
# 19571_2   28            429            69.02473                533  28.0
# 14 13.500000      xstream

# RESULTS DE NOVO:
# T-REKS:
# ["QALAI-","QAV-L-","QAA-Q-","QAV-MG"]
# ["TALE--ALLCPSATVTQ-VCSNPPCETHETGTTN","TATT--SN-AGSAQR--V-CSNPPCETHETGTTH","TATTATSN-GGTGQP--EGGQQPPA-----GRPC"]
# ["TSTVGQQNGSVVRVC-SNPPCETHETGTTNTATTA","TSNMAGQHG-----C-SNPPCETHETGTTNTATTA","MSSVGANHQ-----RDARRACAAG-TPAVIRISVA"]
# TRUST
# ["T--AATA-TS-------PTP-NPvP-SVpANP-----P---KSPAP--A---A---------AA--PAVQPLTQVG----IT-L--L-PQAAP-AP-------PT--TT-TIQVLP","T--VPGS-SI-------SVP-TA-A-RT-QGV-----P--AVLKVT--GpqaT---------TG--TPLVTMRPAS----QAGK--ApVTVTS-LPagvrmvvPT--QSAQ-GTV-","T--LITA-PS-------GVE-AQ-P--V-HDL-----P--VSILAS--P---T---------TEqpTATVTIADSGqgdvQPGT--V-TLVCS-NP-------PC--ETHETGTTN","T--ATTTvVA-------NLG-GH-P-QP-TQV-----QfvCDRQEA--A------------------ASLVTSTVGq---QNGS--V-VRVCS-NP-------PC--ETHETGTTN","T--ATTA-TS-------NMA-GQ-H-GC-SNP-----P--CETHET--G---T---------TN--TATTAMSSVG----ANHQrdA-RRACA------------------AGTPA","VirISVA-TG-------ALE-AA-Q-G--SKS-----Q--CQTRQT--S---A---------TS--TTMTVMAT-G----AP---------CSaGPllg----PS--MAREPGGRS","P--AFVQ-LA-------PLS-SK-V-RL-SSPsikdlP--AGRH-S--H---A---------VS--TAAMTRSSVG-----AGE--P-RMA----P-------VC--ESLQGGSPS","T--TVTV-TAleallcpSAT-VT-Q-VC-SNP-----P--CETHET--G---T---------TN--TATT-----S----NAGS--A-QRVCS-NP-------PC--ETHETGTTH","T--ATTA-TS-------NGGtGQ-PeGG-QQPpagr-P--CETHQ-------T---------TS--TGTTMSVSV-------GA--L-L------P----------------DATS","S--HRTV-ES-------GLE-------V-AAA-----P--SVTPQA--G---T--------------------ALL----APFP--T-QRVCS-NP-------PC--ETHETGTTH","R--AVTT-VT-------QST-PV-P-GP-SVP-----P--PEELQVspG---PrqqlpprqlLQ--SASTAL--MG----ESAE--V-LSASQ-TPel-----PAavDLSSTGEPS"]
# ["SN-PATrmLKTA-AAQVGTSVSSATNTSTRPIITVHK----SGTVTVAQQAQVVTTV----VGG--VTKTITLV","KS-PIS--VPGG-SALISNLGKVMSVVQTKPVQTSAV----TGQASTGPVTQIIQTKGPLPAGT--ILKLVTSA","DGkPTT--IIT--TTQASGAGTKPTILGISSVSPST-----TKPGTTTIIKTIPMSAIITQAGAtgVTSS-PGI","KS-PIT--IITT-KVMTSGTGAPAKIITAVP---KIA----TGHGQQG-VTQVVLKGAPGQPGT--ILRTV---","---PM-----GG-VRLVT----PVTVSAVKPAVTTLVvkgtTGVTTLGTVTGTVSTSLAGAGGH--STS--ASL","AT-PIT--TLGTiATLSSQVINPTAITVSA-AQTTLT----AAGGLTTPTITMQPVSQPTQ-------------"]
# ["WKRVVG----W-SGPVPRPRHGHRAVA------IKELIVVFGggneG--------------I-------------------VDELhvYNTA----TNQWF","IPAVRGDIP------PGCAAYG--FVC------DGTRLLVFG----Gmveygk--------Y-------------------SNDL--YELQ----ASRWE","WKRLKAKTP-K-NGPPPCPRLGHSFSL------VGNKCYLFG----GlandsedpknniprY-------------------LNDL--YILElrpgSGVVA","W-----DIPiT-YGVLPPPRESHTAVVytekdnKKSKLVIYG----Gmsgc----------R-------------------LGDL--WTLD----IDTLT","WNK-----P-SlSGVAPLPRSLHSATT------IGNKMYVFG----G--------------WvplvmddvkvathekewkcTNTL--ACLN----LDTMA","WETILMDTL-E-DNI-PRARAGHCAVA------INTRLYIWS----G-----------------------------------------------------"]
# XSTREAM
# ["CSNPPCETHETGTTNTATTATSNMAG-QHG","CSNPPCETHETGTTNTATTAMS-SVGANHQ","CSNPPCETHETGTTNTATTAMSMAGNHG","CSNPPCETHETGTTNTATTAMS-MAG-NHG"]
# ["VCSNPPCETHETGTTNTA-T-TSNAGSAQR","VCSNPPCETHETGTTHTATTATSNGGTGQP","VCSNPPCETHETGTTHTATTSNAGSAQP","VCSNPPCETHETGTTHTA-T-TSNAGSAQP"]
# ["GTVT","GTVS","GTVS","GTVS"]
# ["AAAE","AAAQ","AAA","AAAE","AAAE"]
# ["ATA","ATA","ATA","ATA"]
# ["EGQ","EGQ","EGQ","EGQ"]
# ["A","A","A","A","A","A","A","A"]


TEST_FAA_FILE = "HIV-1_388796.faa"
TEST_HMM_FILE = "zf-CCHC.hmm"
TEST_DENOVO_PARAMETERS = {"detection": {"detectors": ["XSTREAM", "T-REKS"]}}
TEST_SCORE = "phylo_gap01"
TEST_RESULT_SEQ1 = [3,
                    [["LFNSTKLE", "LFNSST-N"],
                     ["GDII", "GDIR"], ["FLG", "FLG"]],
                    3,
                    [["FNCGG-EF", "FYCNTSNL", "FNSTKLEL", "FNSST-NL"],
                        ["GDII", "GDIR"], ["FLG", "FLG"]]]


@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.mark.slow
def test_MBE_2014_workflow(path):
    # The Schaper et al. (MBE, 2014) workflow is tested on a single sequence.

    test_lSeq = sequence.Sequence.create(
        os.path.join(
            path,
            TEST_FAA_FILE_MBE_2014),
        input_format='fasta')
    test_seq = test_lSeq[0]

    # Information on sequence domains (here: Pfam) in this sequence are added.
    test_pfam_hmm = [
        hmm.HMM.create(
            input_format='hmmer',
            file=os.path.join(
                path,
                i)) for i in TEST_HMM_FILES_MBE_2014]

    # The sequence is searched for tandem repetitions of the Pfam domain in
    # the sequence
    test_pfam_list = test_seq.detect(lHMM=test_pfam_hmm)
    assert len(test_pfam_list.repeats) == 2

    # Pfam TRs with n_effective < 3.5 are discarded.
    test_pfam_list = test_pfam_list.filter("attribute", "n_effective", "min", 3.5)
    assert len(test_pfam_list.repeats) == 2

    # de novo detection methods (Trust, T-reks, Xstream, HHrepID) are used to
    # search the
    test_denovo_list = test_seq.detect(denovo=True, **TEST_DENOVO_PARAMETERS)
    # When Trust is part of the detectors, the number of found repeats may
    # differ between runs...
    assert len(test_denovo_list.repeats) == 10

    # De novo TRs with dTR_units (divergence) > 0.8; n_effective < 2.5; l < 10
    # or pvalue "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions" >0.01
    # are discarded.
    test_denovo_list = test_denovo_list.filter(
        "pvalue",
        TEST_SCORE_MBE_2014,
        0.01)
    assert len(test_denovo_list.repeats) == 10
    test_denovo_list = test_denovo_list.filter(
        "divergence",
        TEST_SCORE_MBE_2014,
        0.8)
    assert len(test_denovo_list.repeats) == 10
    test_denovo_list = test_denovo_list.filter("attribute", "n_effective", "min", 2.5)
    assert len(test_denovo_list.repeats) == 5
    test_denovo_list = test_denovo_list.filter("attribute", "l_effective", "min", 10)
    assert len(test_denovo_list.repeats) == 2

    # De novo TRs were remastered with HMM
    test_denovo_hmm = [
        hmm.HMM.create(
            input_format='repeat',
            repeat=iTR) for iTR in test_denovo_list.repeats]
    test_denovo_list_remastered = test_seq.detect(lHMM=test_denovo_hmm)
    assert len(test_denovo_list_remastered.repeats) == 2

    # pvalue "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions" > 0.1
    # are discarded.
    test_denovo_list_remastered = test_denovo_list_remastered.filter(
        "pvalue",
        TEST_SCORE_MBE_2014,
        0.1)

    # De novo TRs were filtered (n_effective < 3.5 are discarded.)
    test_denovo_list_remastered = test_denovo_list_remastered.filter(
        "attribute",
        "n_effective",
        "min",
        3.5)
    assert len(test_denovo_list_remastered.repeats) == 2

    # De novo TRs overlapping with a Pfam TR were filtered
    test_denovo_list_remastered = test_denovo_list_remastered.filter(
        "none_overlapping_fixed_repeats",
        test_pfam_list,
        "shared_char")
    assert len(test_denovo_list_remastered.repeats) == 2

    # Remaining De novo TRs were clustered for overlap (common ancestry). Only best =
    # lowest p-Value and lowest divergence were retained.
    test_denovo_list_remastered = test_denovo_list_remastered.filter(
        "none_overlapping", ("common_ancestry", None), [
            ("pvalue", TEST_SCORE_MBE_2014), ("divergence", TEST_SCORE_MBE_2014)])
    assert len(test_denovo_list_remastered.repeats) == 1

    # Merge remaining set of de novo and Pfam TRs.
    test_entire_set = test_pfam_list + test_denovo_list_remastered
    assert len(test_entire_set.repeats) == 3

    # Write result set of Pfam TRs
    # test_entire_set.write(format = "tsv, ...")
