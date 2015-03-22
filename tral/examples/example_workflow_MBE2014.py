
'''
Implementation of the workflow used in :
Schaper,E. et al. (2014) Deep conservation of human protein tandem repeats within the eukaryotes. Molecular Biology and Evolution. 31, 1132–1148 .
'''

import logging
import logging.config
import os

from tral.paths import *
from tral.paths import PACKAGE_DIRECTORY

from tral.sequence import repeat_detection_run, sequence
from tral.hmm import hmm

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')

TEST_FAA_FILE_MBE_2014 = os.path.join(PACKAGE_DIRECTORY,"test","P51610.fasta")
TEST_HMM_FILES_MBE_2014 = [os.path.join(PACKAGE_DIRECTORY,"test","Kelch_1.hmm"), os.path.join(PACKAGE_DIRECTORY,"test","Kelch_2.hmm")]
TEST_SCORE_MBE_2014 = "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions"

def path():
    """Return the path to the test data files.
    """
    return os.path.join(os.path.abspath('.'), 'tandemrepeats', 'test')

def sample_MBE_2014_pipeline():
    # The Schaper et al. (MBE, 2014) pipeline is tested on a single sequence.

    test_lSeq = sequence.Sequence.create(os.path.join(path(), TEST_FAA_FILE_MBE_2014), format = "fasta")
    test_seq = test_lSeq[0]

    # Information on sequence domains (here: Pfam) in this sequence are added.
    test_pfam_hmm = [hmm.HMM.create( format = "hmmer", file = os.path.join(path(),i) ) for i in TEST_HMM_FILES_MBE_2014]

    # The sequence is searched for tandem repetitions of the Pfam domain in the sequence
    test_pfam_list = test_seq.detect(lHMM = test_pfam_hmm)
    assert len(test_pfam_list.repeats) == 2

    # Pfam TRs with nD < 3.5 are discarded.
    test_pfam_list = test_pfam_list.filter("attribute", "nD", "min", 3.5)
    assert len(test_pfam_list.repeats) == 2

    # de novo detection methods (Trust, T-reks, Xstream, HHrepID) are used to search the
    # INSERT OWN PARAMTERS USING: test_denovo_list = test_seq.detect(denovo = True, **TEST_DENOVO_PARAMETERS)
    test_denovo_list = test_seq.detect(denovo = True)
    # When Trust is part of the detectors, the number of found repeats may differ between runs...
    assert len(test_denovo_list.repeats) == 10

    # De novo TRs with dTR_units (divergence) > 0.8; nD < 2.5; l < 10 or
    # pValue "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions" > 0.01 are discarded.
    test_denovo_list = test_denovo_list.filter("pValue", TEST_SCORE_MBE_2014, 0.01)
    assert len(test_denovo_list.repeats) == 10
    test_denovo_list = test_denovo_list.filter("divergence", TEST_SCORE_MBE_2014, 0.8)
    assert len(test_denovo_list.repeats) == 10
    test_denovo_list = test_denovo_list.filter("attribute", "nD", "min", 2.5)
    assert len(test_denovo_list.repeats) == 5
    test_denovo_list = test_denovo_list.filter("attribute", "l", "min", 10)
    assert len(test_denovo_list.repeats) == 2

    # De novo TRs were remastered with HMM
    test_denovo_hmm = [hmm.HMM.create(repeat = iTR) for iTR in test_denovo_list.repeats]
    test_denovo_list_remastered = test_seq.detect(lHMM = test_denovo_hmm)
    assert len(test_denovo_list_remastered.repeats) == 2

    # pValue "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions" > 0.1 are discarded.
    test_denovo_list_remastered = test_denovo_list_remastered.filter("pValue", TEST_SCORE_MBE_2014, 0.1)

    # De novo TRs were filtered (nD < 3.5 are discarded.)
    test_denovo_list_remastered = test_denovo_list_remastered.filter("attribute", "nD", "min", 3.5)
    assert len(test_denovo_list_remastered.repeats) == 2

    # De novo TRs overlapping with a Pfam TR were filtered
    test_denovo_list_remastered = test_denovo_list_remastered.filter("none_overlapping_fixed_repeats", test_pfam_list, "shared_char")
    assert len(test_denovo_list_remastered.repeats) == 2

    # Remaining De novo TRs were clustered for overlap (common ancestry). Only best =
    # lowest p-Value and lowest divergence were retained.
    test_denovo_list_remastered = test_denovo_list_remastered.filter("none_overlapping", ["common_ancestry"], {"pValue":TEST_SCORE_MBE_2014, "divergence":TEST_SCORE_MBE_2014})
    assert len(test_denovo_list_remastered.repeats) == 1

    # Merge remaining set of de novo and Pfam TRs.
    test_entire_set  = test_pfam_list + test_denovo_list_remastered
    assert len(test_entire_set.repeats) == 3

    # Write result set of Pfam TRs
    #test_entire_set.write(format = "tsv,...")


if __name__ == "__main__":

    sample_MBE_2014_pipeline()
