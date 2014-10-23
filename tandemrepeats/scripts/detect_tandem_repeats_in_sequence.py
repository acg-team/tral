#!python

import argparse
import itertools
import json
import logging
import logging.config
import os
import pickle
import sys

from tandemrepeats import configuration
from tandemrepeats.paths import *
from tandemrepeats.repeat_list import repeat_list
from tandemrepeats.sequence import sequence
from tandemrepeats.hmm import hmm

logging.config.fileConfig(os.path.join(CODEROOT,'tandemrepeats','data','logging.ini'))
log = logging.getLogger('root')

c = configuration.Configuration.Instance()
config = c.config

# Shift REPEAT_LIST_TAG to configuration file and rename?
REPEAT_LIST_TAG = "all"
DE_NOVO_ALL_TAG = "denovo_all"
PFAM_ALL_TAG = "pfam_all"
DE_NOVO_TAG = "denovo"
PFAM_TAG = "pfam"
DE_NOVO_REFINED_TAG = "denovo_refined"
DE_NOVO_FINAL_TAG = "denovo_final"

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
    lHMM = [hmm_ID for iS in lSequence for hmm_ID in iS.get_annotation('PFAM')]
    infoNRuns = len(lHMM)
    log.debug("{} Viterbi runs need to be performed.".format(infoNRuns))
    lHMM = set(lHMM)
    infoNHMM = len(lHMM)
    log.debug("These derive from {} independent HMMs.".format(infoNHMM))
    dHMM = {hmm_ID: hmm.HMM.create(format = "pickle", file = os.path.join(hmm_dir, hmm_ID + ".pickle"))
                for hmm_ID in lHMM}

    dTR = {}
    for iS in lSequence:
        log.debug("".format(iS.get_annotation('PFAM')))
        if iS.get_annotation('PFAM'):
            iRL = iS.detect([dHMM[hmm_ID] for hmm_ID in iS.get_annotation('PFAM')])
            for iTR, hmm_ID in zip(iRL.repeats, iS.get_annotation('PFAM')):
                iTR.model = hmm_ID
            dTR[iS.id] = iRL
        else:
           dTR[iS.id] = None

    with open(result_file, 'wb') as fh:
        pickle.dump(dTR, fh)

    print("\n" + "\t".join(["infoNRuns", "infoNHMM"]))
    print("\n" + "\t".join([infoNRuns, infoNHMM]))
    print("DONE")


def annotate_de_novo(sequences_file, result_file, detector = None):
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
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    if detector:
        detection_parameters = {"detection": {"lFinders": [detector]}}
    else:
        detection_parameters = {}
    log.debug("detection_parameters: {}".format(detection_parameters))

    dRL = {}
    for iS in lSequence:
        iRL = iS.detect(denovo = True, **detection_parameters)
        log.debug(iRL.repeats)
        for iTR in iRL.repeats:
            iTR.TRD = detector
        dRL[iS.id] = iRL

    with open(result_file, 'wb') as fh:
        pickle.dump(dRL, fh)

    print("DONE")


def calculate_significance(repeat_file, result_file, **kwargs):
    ''' Calculate the statistical significance for all repeats in ``repeat_file``.

     Calculate the statistical significance for all repeats in ``repeat_file``.

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
        if iRL:
            for iTR in iRL.repeats:
                iTR.calculate_pValues(**kwargs)

    with open(result_file, 'wb') as fh:
        pickle.dump(dRL, fh)

    print("DONE")


def merge_and_basic_filter(sequences_file, repeat_files, result_file, **kwargs):

    ''' Merge TR annotations from several sources and perform basic filtering.

    Merge TR annotations from several sources and perform basic filtering based on the
    number of repeat units and the statistical significance of the tandem repeats.

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

    ..ToDo: Input BASIC_FILTER as PARAMETERS
    '''

    basic_filter = config['filter']['basic']
    basic_filter_tag = config['filter']['basic']['tag']

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

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
        iS.set_repeat_list(dRL_all[iS.id], REPEAT_LIST_TAG)
        if iS.dRepeat_list[REPEAT_LIST_TAG]:
            denovo_repeat_list = repeat_list.Repeat_list([i for i in iS.dRepeat_list[REPEAT_LIST_TAG].repeats if hasattr(i, "TRD")])
            pfam_repeat_list = repeat_list.Repeat_list([i for i in iS.dRepeat_list[REPEAT_LIST_TAG].repeats if not hasattr(i, "TRD")])
        else:
            denovo_repeat_list = iS.dRepeat_list[REPEAT_LIST_TAG]
            pfam_repeat_list = iS.dRepeat_list[REPEAT_LIST_TAG]
        iS.set_repeat_list(denovo_repeat_list, DE_NOVO_ALL_TAG)
        iS.set_repeat_list(pfam_repeat_list, PFAM_ALL_TAG)

    for iS in lSequence:
        rl_tmp = iS.dRepeat_list[REPEAT_LIST_TAG]
        if iS.dRepeat_list[REPEAT_LIST_TAG]:
            for iB in basic_filter['dict'].values():
                rl_tmp = rl_tmp.filter(**iB)
        else:
            rl_tmp = iS.dRepeat_list[REPEAT_LIST_TAG]
        iS.set_repeat_list(rl_tmp, basic_filter_tag)
        iS.set_repeat_list(rl_tmp.intersection(iS.dRepeat_list[PFAM_ALL_TAG]), PFAM_TAG)
        iS.set_repeat_list(rl_tmp.intersection(iS.dRepeat_list[DE_NOVO_ALL_TAG]), DE_NOVO_TAG)

    with open(result_file, 'wb') as fh:
        pickle.dump(lSequence, fh)

    print("DONE")


def calculate_overlap(sequences_file, result_file, lOverlap_type, **kwargs):
    ''' Calculate the overlap of TR annotations from several sources.

    Calculate the overlap of TR annotations from several sources.

    Args:
        sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
        result_file (str): Path to the result file.
        kwargs (dict): A dictionary of parameters for the applied significance test.


     Raises:
        Exception: If any of the pickles in ``sequences_file`` cannot be loaded
    '''

    basic_filter_tag = config['filter']['basic']['tag']

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    for iS,iO in itertools.product(lSequence, lOverlap_type):
        iS.dRepeat_list[basic_filter_tag].cluster(overlap_type = iO)

    for iS in lSequence:
        # Perform common ancestry overlap filter and keep PFAMs
        iC = {"func_name": "none_overlapping_fixed_repeats", "rl_fixed": iS.dRepeat_list[PFAM_TAG], "overlap_type": "common_ancestry"}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**iC)

        # Choose only the most convincing de novo TRs
        iC = {"func_name": "none_overlapping", "overlap": ("common_ancestry", None), "lCriterion": [("pValue", "phylo_gap01"), ("divergence", "phylo_gap01")]}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**iC)

    with open(result_file, 'wb') as fh:
        pickle.dump(lSequence, fh)

    print("DONE")


def refine_denovo(sequences_file, result_file):
    ''' Refine denovo TRs.

     Refine denovo TRs. Calculate significance. Check location in sequence. Discard insignificant.

     Args:
         sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
         result_file (str): Path to the result file.

     Raises:
        Exception: If the pickle ``sequences_file`` cannot be loaded
    '''

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    for iS in lSequence:
        log.debug(iS.id)
        denovo_final = []
        denovo_refined = []
        for iTR in iS.dRepeat_list[DE_NOVO_TAG].repeats:
            # Create HMM from TR
            denovo_hmm = hmm.HMM.create(format = 'repeat', repeat = iTR)
            # Run HMM on sequence
            denovo_refined_rl = iS.detect(lHMM = [denovo_hmm])
            if denovo_refined_rl:
                iTR_refined = denovo_refined_rl.repeats[0]
                iTR_refined.TRD = iTR.TRD + "_refined"
                denovo_refined.append(iTR_refined)
                # Check whether new and old TR overlap. Check whether new TR is significant. If not both, put unrefined TR into final.
                if not repeat_list.two_repeats_overlap("shared_char", iTR, iTR_refined) and not iTR_refined.pValue("phylo_gap01") < 0.1:
                    denovo_final.append(iTR)
                else:
                    denovo_final.append(iTR_refined)
        iS.set_repeat_list(repeat_list.Repeat_list(denovo_refined), DE_NOVO_REFINED_TAG)
        iS.set_repeat_list(repeat_list.Repeat_list(denovo_final), DE_NOVO_FINAL_TAG)

    with open(result_file, 'wb') as fh:
        pickle.dump(lSequence, fh)

    print("DONE")


def serialize_annotations(sequences_dir, result_file, format):
    ''' Serialize tandem repeat annotation results

     Serialize tandem repeat annotation results

     Args:
         sequences_dir (str): Path to the dir containing pickle file of lists of ``Sequence``
            instances.
         result_file (str): Path to the result file.
         format (str): Serialization format name, e.g. tsv.

     Raises:
        Exception: If the pickles in ``sequences_dir`` cannot be loaded

    ..ToDo: Serialisation format currently not used..
    '''

    lFiles = [file for file in os.listdir(sequences_dir) if file.endswith(".pickle")]

    with open(result_file, 'w') as fh:
        for iFile in Files:
            sequences_file = os.path.join(sequences_dir, iFile)
            try:
                with open(sequences_file, 'rb') as fh:
                    lSequence = pickle.load(fh)
            except:
                raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

            for result_file in lSequence:
                for iTR in iS.dRepeat_list[PFAM_TAG].repeats:
                    data = [str(i) for i in [iS.id, " ".join(iTR.msa), iTR.begin, iTR.pValue("phylo_gap01"), iTR.lD, iTR.n, iTR.nD, iTR.model]]
                    result_file.write("\t".join(data))
                for iTR in iS.dRepeat_list[DE_NOVO_TAG].repeats:
                    data = [str(i) for i in [iS.id, " ".join(iTR.msa), iTR.begin, iTR.pValue("phylo_gap01"), iTR.lD, iTR.n, iTR.nD, iTR.TRD]]
                    result_file.write("\t".join(data))

    print("DONE")


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

    # method = getattr(sys.modules[__name__], pars["method_name"])
    # method()
    if pars["method"] == "annotate_de_novo":
        annotate_de_novo(pars["input"], pars["output"], detector = pars['detectors'])
    elif pars["method"] == "annotate_TRs_from_hmmer":
        annotate_TRs_from_hmmer(pars["input"], pars["hmm"], pars["output"])
    elif pars["method"] == "calculate_significance":
        calculate_significance(pars["input"], pars["output"])
    elif pars["method"] == "merge_and_basic_filter":
        merge_and_basic_filter(pars["input"], pars['repeat_files'], pars["output"])
    elif pars["method"] == "calculate_overlap":
        calculate_overlap(pars["input"], pars["output"], pars["overlap_type"])
    elif pars["method"] == "refine_denovo":
        refine_denovo(pars["input"], pars["output"])
    elif pars["method"] == "serialize_annotations":
        serialize_annotations(pars["input"], pars["output"], )

def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('input', metavar='input_file', type=str,
                       help='The path to the input file.')
    parser.add_argument('output', metavar='output_file', type=str,
                       help='The path to the output file.')
    parser.add_argument('method', metavar='method_name', type=str,
                       help='The name of the method to be executed.')
    parser.add_argument('-seq','--sequence_type', type=str,
                       help='The sequence type: -seq AA or -seq DNA')
    parser.add_argument('-rep','--repeat_files', nargs='+', type=str,
                       help='The repeat files')
    parser.add_argument('-f','--format', type=str,
                       help='The serialization format, e.g. tsv')
    parser.add_argument('-ov','--overlap_type', nargs='+', type=str,
                       help='The overlap type, e.g. "common_ancestry shared_char"')
    parser.add_argument('-d', '--detectors', nargs='+', type=str,
                        help='The de novo tandem repeat detectors. For example: -d T-REKS XSTREAM')
    parser.add_argument('-hmm', '--hmm', type=str,
                        help='The path to the dir with HMMER pickles')
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