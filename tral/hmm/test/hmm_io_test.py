import os
import pytest

from tral.hmm.hmm import HMM
from tral.hmm import hmm_io
from tral.repeat.repeat import Repeat


NUM_STATES = 8
AMINOACIDS = list('ACDEFGHIKLMNPQRSTVWY')
NUM_AMINOACIDS = len(AMINOACIDS)
NUM_TRANSITIONS = 7

CARCINUSTATIN_ID = 'PF08261.7'
WRONG_CARCINUSTATIN_ID = 'PF05932'
SHORT_CARCINUSTATIN_ID = 'PF08261'

# HMM data dictionary key names
INSERTION_EMISSIONS = 'insertion_emissions'
EMISSIONS = 'emissions'
TRANSITIONS = 'transition'
COMPO_KEY_NAME = 'COMPO'
LETTER_KEY_NAME = 'letters'
ID_KEY_NAME = 'id'

# Test file names
TEST_FILE_WITH_ID = 'carcinustatin.hmm'
TEST_FILE_WITHOUT_ID = 'carcinustatin_no_id.hmm'


TEST_REPEAT_MSA_SINGLE = ["A","A","A"]


@pytest.fixture
def path():
    """Return the path to the test data files.
    """
    return os.path.dirname(os.path.abspath(__file__))

@pytest.mark.no_external_software_required
def test_single_hmm_with_id_read(path):
    test_dict_list = list(HMM.read(os.path.join(path, TEST_FILE_WITH_ID)))
    assert len(test_dict_list) == 1
    test_dict = test_dict_list[0]

    compare_carcinustatin(test_dict)

    assert test_dict[ID_KEY_NAME] == CARCINUSTATIN_ID

@pytest.mark.no_external_software_required
def test_single_hmm_without_id_read(path):
    test_dict_list = list(HMM.read(os.path.join(path, TEST_FILE_WITHOUT_ID)))
    assert len(test_dict_list) == 1
    test_dict = test_dict_list[0]

    compare_carcinustatin(test_dict)

    assert test_dict[ID_KEY_NAME] is None


@pytest.mark.no_external_software_required
def test_single_hmm_no_id_with_query_read(path):
    test_dict_list = list(HMM.read(os.path.join(path, TEST_FILE_WITHOUT_ID),
                              id=CARCINUSTATIN_ID))
    assert len(test_dict_list) == 0


@pytest.mark.no_external_software_required
def test_single_hmm_with_wrong_query_read(path):
    test_dict_list = list(HMM.read(os.path.join(path, TEST_FILE_WITH_ID),
                              id=WRONG_CARCINUSTATIN_ID))
    assert len(test_dict_list) == 0


@pytest.mark.no_external_software_required
def test_single_hmm_with_short_query_read(path):
    test_dict_list = list(hmm_io.read(os.path.join(path, TEST_FILE_WITH_ID),
                              id=SHORT_CARCINUSTATIN_ID))
    assert len(test_dict_list) == 1
    test_dict = test_dict_list[0]

    compare_carcinustatin(test_dict)

    assert test_dict[ID_KEY_NAME] == CARCINUSTATIN_ID

@pytest.mark.no_external_software_required
def compare_carcinustatin(test_dict):
    for i in range(0, NUM_STATES):
        assert str(i + 1) in test_dict.keys()
        assert INSERTION_EMISSIONS in test_dict[str(i + 1)].keys()
        assert EMISSIONS in test_dict[str(i + 1)].keys()
        assert TRANSITIONS in test_dict[str(i + 1)].keys()
        assert len(test_dict[str(i + 1)][INSERTION_EMISSIONS]) ==\
            NUM_AMINOACIDS
        assert len(test_dict[str(i + 1)][EMISSIONS]) ==\
            NUM_AMINOACIDS
        assert len(test_dict[str(i + 1)][TRANSITIONS]) == NUM_TRANSITIONS

    assert INSERTION_EMISSIONS in test_dict[COMPO_KEY_NAME].keys()
    assert EMISSIONS in test_dict[COMPO_KEY_NAME].keys()
    assert TRANSITIONS in test_dict[COMPO_KEY_NAME].keys()
    assert len(test_dict[COMPO_KEY_NAME][INSERTION_EMISSIONS]) == NUM_AMINOACIDS
    assert len(test_dict[COMPO_KEY_NAME][EMISSIONS]) == NUM_AMINOACIDS
    assert len(test_dict[COMPO_KEY_NAME][TRANSITIONS]) == NUM_TRANSITIONS

    assert len(set(test_dict[LETTER_KEY_NAME]) ^ set(AMINOACIDS)) == 0

    assert ID_KEY_NAME in test_dict.keys()
