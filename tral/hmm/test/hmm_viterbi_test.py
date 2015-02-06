import os
import pytest

from tandemrepeats.hmm.hmm import HMM
from tandemrepeats.hmm.hmm_viterbi import *
from tandemrepeats.repeat import repeat
from tandemrepeats.repeat import repeat_align


notfixed = pytest.mark.notfixed

def test_hmm_path_to_non_aligned_tandem_repeat_units():

    TEST = {"Double": ["AAAA", ["M1","M2","M1","M2"], 2, ["AA","AA"]],
            "Double_converted": ["ABCD", ["M2","M1","M2","M1"], 2, ["AB","CD"]],
            "Double_converted_complex": ["NNAAAACC", ["N","N","M2","M1","M2","M1","C","C"], 2, ["AA","AA"]],
            "Double_converted_complex_insertions": ["NNAIAAACC", ["N","N","M2","I2","M1","M2","M1","C","C"], 2, ["AIA","AA"]],
            "Single": ["AAAA", ["M1","M1","M1","M1"], 1, ["A","A","A","A"]],
            "Single_Complex": ["NNAAIIIAA", ["N","N","M1","M1","I1","I1","I1","M1","M1"], 1, ["A","AIII","A","A"]],
            "Long": ["GYRADKLADKLADKL", ["N","N","N","M1","M2","M3","M4","M1","M2","M3","M4","M1","M2","M3","M4"], 4, ["ADKL","ADKL","ADKL"]],
            }
    for test, p in TEST.items():
        test_repeat_msa = hmm_path_to_non_aligned_tandem_repeat_units(sequence = p[0], path = p[1], lD = p[2])
        assert test_repeat_msa == p[3]

def test_viterbi():

    # {Test_name: [Original_TR_MSA, Sequence, Viterbi_path, Refined_TR_MSA], ... }
    TEST = {"Single": [["A","A","A"], "AAAAAA", ["M1","M1","M1","M1","M1","M1"], ["A","A","A","A","A","A"]],
        "Double": [["AA","AA"], "AAAAAA", ["M1","M2","M1","M2","M1","M2"], ["AA","AA","AA"]],
        "Long": [["ADKL","ADKL"], "GYRADKLADKLADKL", ["N","N","N","M1","M2","M3","M4","M1","M2","M3","M4","M1","M2","M3","M4"], ["ADKL","ADKL","ADKL"]]
        }

    for test, p in TEST.items():
        test_repeat = repeat.Repeat(msa = p[0])
        test_hmm = HMM.create(format = "repeat", repeat = test_repeat)

        for iHMM in [test_hmm]:
            # Detect TRs on self.seq with hmm using the Viterbi algorithm.
            most_likely_path = iHMM.viterbi(p[1])
            assert type(most_likely_path) == list
            assert most_likely_path == p[2]

            unaligned_msa = hmm_path_to_non_aligned_tandem_repeat_units(p[1], most_likely_path, iHMM.lD)
            assert unaligned_msa == p[3]

            aligned_msa = repeat_align.realign_repeat(unaligned_msa)
            assert aligned_msa == p[3]

# The following test functions were copy-pasted from hmm_viterbi.py. They need to be either
# strongly adapted, or - if not useful - erased.
@notfixed
def test_conversion_multiple():
    """ Test <hmm_path_to_maximal_complete_tandem_repeat_units>.
        However, no positive set is defined right now.
    """
    assert 0, "Test not fixed"
    lD = 3
    lPaths = [['N','N','N', 'M0','M1','M2','M0','I1','M1','M2','M0','M1','M2','N','N','N'],
            ['N','N','N','M1','M2','M0','M1','M2','M0','M1','M2', 'M0', 'N','N','N'],
            ['N','N','N','M2','M0','M1','M2','M0','M1','M2','M0','M1', 'N','N','N']]
    lSequence = ['XXXABCADBCABCXXX','XXXBCABCABCAXXX','XXXCABCABCABXXX']
    lMSA = hmm_path_to_maximal_complete_tandem_repeat_units(lSequence, lPaths, 3)

@notfixed
def test_create_viterbi():
    """ This test needs to be fixed, the HMM cannot be initialised in this way.
    """
    assert 0, "Test not fixed"
    my_hmm = HMM()
    my_hmm.states = ['H', 'F']

    # Initialisation
    my_hmm.p_0 = {"H": 0.6, "F": 0.4}

    # Feed Values to p_t
    #my_hmm.p_t["START"] = {"H": 0.6, "F": 0.4}
    my_hmm.p_t["H"] = {"H": 0.7, "F": 0.3}
    my_hmm.p_t["F"] = {"H": 0.4, "F": 0.6}

    # emissions
    my_hmm.emissions = ["D", "C", "N"]

    # emission probabilities
    my_hmm.p_e['H'] = {"D":0.1, "C":0.4, "N":0.5}
    my_hmm.p_e['F'] = {"D":0.6, "C":0.3, "N":0.1}

    my_viterbi = Viterbi(my_hmm, "NCD")
    if my_viterbi.viterbi() == "HHF":
        print("Test SUCCESSFUL")
    else:
        print("Test FAILED")

@notfixed
def test_init_with_repeat():
    """ This test needs to be fixed.
    """
    assert 0, "Test not fixed"
    my_TR = Repeat(begin = 0, msa = ['A-G', 'ACG', 'ACG'], sequence_type = 'DNA')
    my_hmm = HMM.create(repeat=my_TR)
    from . import sequence
    my_sequence = sequence.Sequence()
    my_viterbi = Viterbi(my_hmm, my_sequence.sequence)
    v = my_viterbi.viterbi()
    "".join(v)
