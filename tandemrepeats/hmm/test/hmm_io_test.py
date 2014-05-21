import logging
import os

from tandemrepeats.hmm.hmm import HMM

NUM_STATES = 8
AMINOACIDS = list('ACDEFGHIKLMNPQRSTVWY')
NUM_AMINOACIDS = len(AMINOACIDS)

# HMM data dictionary key names
INSERTION_EMISSIONS = 'insertion_emissions'
EMISSIONS = 'emissions'
TRANSITIONS = 'transition'
COMPO_KEY_NAME = 'COMPO'
LETTER_KEY_NAME = 'letters'

def test_single_hmm_read():

    test_dict_list = HMM.read(os.path.join(os.path.abspath('.'), 'hmm', 'test',
                                           'carcinustatin.hmm'))
    assert len(test_dict_list) == 1
    test_dict = test_dict_list[0]

    for i in range(0, NUM_STATES):
        assert str(i + 1) in test_dict.keys()
        assert INSERTION_EMISSIONS in test_dict[str(i + 1)].keys()
        assert EMISSIONS in test_dict[str(i + 1)].keys()
        assert TRANSITIONS in test_dict[str(i + 1)].keys()
        assert len(test_dict[str(i + 1)][INSERTION_EMISSIONS]) ==\
            NUM_AMINOACIDS
        assert len(test_dict[str(i + 1)][EMISSIONS]) ==\
            NUM_AMINOACIDS
        assert len(test_dict[str(i + 1)][TRANSITIONS]) == NUM_STATES - 1

    assert INSERTION_EMISSIONS in test_dict[COMPO_KEY_NAME].keys()
    assert EMISSIONS in test_dict[COMPO_KEY_NAME].keys()
    assert TRANSITIONS in test_dict[COMPO_KEY_NAME].keys()
    assert len(test_dict[COMPO_KEY_NAME][INSERTION_EMISSIONS]) == NUM_AMINOACIDS
    assert len(test_dict[COMPO_KEY_NAME][EMISSIONS]) == NUM_AMINOACIDS
    assert len(test_dict[COMPO_KEY_NAME][TRANSITIONS]) == NUM_STATES - 1

    assert len(set(test_dict[LETTER_KEY_NAME]) ^ set(AMINOACIDS)) == 0
