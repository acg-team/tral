#!python

import argparse
import json
import logging
import logging.config
import os
import pickle
import sys

from tandemrepeats import configuration
from tandemrepeats.paths import *
from tandemrepeats.sequence import sequence

logging.config.fileConfig(os.path.join(CODEROOT,'tandemrepeats','data','logging.ini'))
log = logging.getLogger('root')

c = configuration.Configuration.Instance()
config = c.config

def annotate_TRs_from_hmmer(sequences_file, hmm_dir, result_file, **kwargs):
    ''' Annotate sequences with TRs from HMMer models.

     Annotate sequences with TRs from HMMer models. Save the annotations in a pickle.

     Args:
         sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
         hmm_dir (str): Path to directory where all HMMs are stored as .pickles
         result_file (str): Path to the result file.

     Raises:
        Exception: If the pickle ``sequences_file`` cannot be loaded
        Exception: if the hmm_dir does not exist
    '''

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    if not os.path.isdir(hmm_dir):
        raise Exception("hmm_dir does not exists: {}".format(hmm_dir))

    # Load all HMM pickles needed for all sequences.
    lHMM = [hmm_ID for iS in lSequence for hmm_ID in iS.lHMM_ID]
    infoNRuns = len(lHMM)
    log.debug("{} Viterbi runs need to be performed.".format(infoNHMM))
    lHMM = set(lHMM)
    infoNHMM = len(lHMM)
    log.debug("These derive from {} independent HMMs.".format(infoNHMM))
    dHMM = {hmm_ID: hmm.HMM.create(format = "pickle", file = os.path.join(hmm_dir, hmm_ID, ".pickle"))
                for iHMM in lHMM}

    dTR = {}
    for iS in lSequence:
        log.debug("".format(iS.lHMM_ID))
        iRL = iS.detect([dHMM[hmm_ID] for hmm_ID in iS.lHMM_ID])
        for iTR, hmm_ID in zip(iRL.repeats, iS.lHMM_ID):
            iTR.model = hmm_ID
        dTR[iS.id] = iRL

    with open(result_file, 'wb') as fh:
        pickle.dump(dTR, fh)

    print("\n" + "\t".join(["infoNRuns", "infoNHMM"]))
    print("\n" + "\t".join([infoNRuns, infoNHMM]))
    print("DONE")


def annotate_de_novo(sequences_file, result_file, detector):
    ''' Annotate sequences with TRs with a de novo TR ``detector``.

     Annotate sequences with TRs with a de novo TR ``detector``.

     Args:
         sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
         detector (str): A tandem repeat de novo detector
         result_file (str): Path to the result file.

     Raises:
        Exception: If the pickle ``sequences_file`` cannot be loaded
    '''

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_File: {}".format(sequences_File))

    detection_parameters = {"detection": {"lFinders": [detector]}}
    log.debug("detection_parameters: {}".format(detection_parameters))

    dRL = {}
    for iS in lSequence:
        iRL = iS.detect(denovo = True, **detection_parameters)
        for iTR in iRL.repeats:
            iTR.TRD = detector
        dRL[iS.id] = iRL

    with open(result_file, 'wb') as fh:
        pickle.dump(dRL, fh)

    print("DONE")


def calculate_significance(repeat_file, result_file, **kwargs):
    ''' Calculate the statistical significance for all repeat in ``repeat_file``.

     Calculate the statistical significance for all repeat in ``repeat_file``.

     Args:
         repeat_file (str): Path to the pickle file containing a dict of ``Repeat``
            instances.
         result_file (str): Path to the result file.
         kwargs (dict): A dictionary of parameters for the applied significance test.

     Raises:
        Exception: If the pickle ``repeat_file`` cannot be loaded
    '''

    try:
        with open(repeat_file, 'rb') as fh:
            dRL = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file repeat_file: {}".format(repeat_file))

    for iRL in dRL.values():
        for iTR in iRL.repeats():
            iTR.calculate_pValues(**kwargs)

    with open(result_file, 'wb') as fh:
        pickle.dump(dRL, fh)

    print("DONE")


def calculate_overlap(sequences_file, repeat_files, result_file, lOverlap_type, **kwargs):
    ''' Calculate the overlap of TR annotations from several sources.

    Calculate the overlap of TR annotations from several sources.

    Args:
        sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
        repeat_files (list of str): Lists of paths to the pickle file containing a dict
            of ``Repeat`` instances.
        result_file (str): Path to the result file.
        kwargs (dict): A dictionary of parameters for the applied significance test.


     Raises:
        Exception: If the pickle ``repeat_file`` cannot be loaded
     Raises:
        Exception: If any of the pickles in ``sequences_file`` cannot be loaded
    '''

    rl_tag = "all"

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_File))

    log.debug("Merging all repeat lists in from ``repeat_files``.")
    dRL_all = {}
    for iRLF in repeat_files:
        try:
            with open(iRLF, 'rb') as fh:
                dRL = pickle.load(fh)
        except:
            raise Exception("Cannot load putative pickle file repeat_list_file: {}".format(iRLF))
        if not dRL_all:
            dRL_all = dRL
        else:
            for iS_ID,iRL in dRL.items():
                dRL_all[iS_ID] += iRL

    log.debug("Append ``repeat_list`` to ``sequence``.")
    for iS in lSequence:
        iS.set_repeat_list(dRL_all[iS.id], rl_tag)

    for iS, iO in zip(lSequence, lOverlap_type):
        iS.dRepeat_list[rl_tag].cluster(overlap_type = iO)

    with open(result_file, 'wb') as fh:
        pickle.dump(lSequence, fh)

    print("DONE")


def filter():
    return True

def refine_denovo():
    return True

def serialize_annotations():
    return True


def main():

    pars = read_commandline_arguments()

    # Update configuration
    if "path" in pars:
        fasta_file = pars["path"]
    if "sequence_type" in pars:
        config["sequence_type"] = pars["sequence_type"]
    if "detectors" in pars:
        config["sequence"]["repeat_detection"][config["sequence_type"]] = pars["detectors"]
    if "significance_test" in pars:
        config["repeat_score"]["score_calibration"] = pars["significance_test"]

    seq = sequence.Sequence.read(fasta_file)[0]
    denovo_list = seq.detect(denovo = True)
    sys.stdout.write(denovo_list.write("tsv") + "\n")

def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('path', metavar='seq.faa', type=str,
                       help='The path to the sequence fasta file.')
    parser.add_argument('-seq','--sequence_type', type=str,
                       help='The sequence type: -seq AA or -seq DNA')
    parser.add_argument('-d', '--detectors', nargs='+', type=str,
                        help='The de novo tandem repeat detectors. For example: -d T-REKS XSTREAM')
    parser.add_argument("-test", "--significance_test", type=str, required=False,
                        help='The significance tests and cut-off values used. For example: -test \'{"phylo_gap01":0.01}\'')
    parser.add_argument("-c", "--cluster", type=str, required=False,
                        help='The overlap definition by which tandem repeats are clustered. For example: -c common_ancestry')
    parser.add_argument("-cf", "--cluster_filter", type=str, required=False,
                        help='The filter by which best representatives are chosen. For example: [TODO]')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()