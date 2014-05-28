import logging
import os
import pytest

from tandemrepeats.hmm.hmm import HMM
from tandemrepeats.hmm.hmm_viterbi import *
from tandemrepeats.repeat.repeat import Repeat

# These test functions were copy-pasted from hmm_viterbi.py. They need to be either
# strongly adapted, or - if not useful - erased.

notfixed = pytest.mark.notfixed

@notfixed
def test_conversion_mutliple():
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
