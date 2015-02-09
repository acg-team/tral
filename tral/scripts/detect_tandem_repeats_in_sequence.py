#!python

import argparse
import datetime
import itertools
import json
import logging
import logging.config
import os
import pickle
import sys

from pyfaidx import Fasta

from tral import configuration
from tral.paths import *
from tral.repeat_list import repeat_list
from tral.sequence import sequence
from tral.hmm import hmm

logging.config.fileConfig(config_file("logging.ini"))
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
FINAL_TAG = "final"

def workflow(sequences_file, hmm_annotation_file, hmm_dir, result_file, max_time, **kwargs):

    ''' Annotate sequences with TRs from multiple sources, test and refine annotations.

     Save the annotations in a pickle.

     Args:
         sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
         hmm_dir (str): Path to directory where all HMMs are stored as .pickles
         result_file (str): Path to the result file.
         max_time (str): Max run time in seconds

     Raises:
        Exception: If the pickle ``sequences_file`` cannot be loaded
        Exception: if the hmm_dir does not exist

     ToDo: Make sure the pValue is automatically calculated!!!!!
    '''

    start = datetime.datetime.now()
    max_time = int(max_time)
    time_interval = 3600
    next_time = 3600

    try:
        lSequence = Fasta(sequences_file)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    if not os.path.isdir(hmm_dir):
        raise Exception("hmm_dir does not exists: {}".format(hmm_dir))

    try:
        with open(hmm_annotation_file, 'rb') as fh:
            dHMM_annotation = pickle.load(fh)
    except:
        raise Exception("Cannot load hmm_annotation_file: {}".format(hmm_annotation_file))

    basic_filter = config['filter']['basic']['dict']
    basic_filter_tag = config['filter']['basic']['tag']

    # Load previous results:
    try:
        with open(result_file, 'rb') as fh:
            dResults = pickle.load(fh)
    except:
        log.debug("Could not load previous results file. Perhaps non existant.")
        dResults = {}


    for iS_pyfaidx in lSequence:

        elapsed_time = (datetime.datetime.now() - start).seconds
        if elapsed_time > max_time or elapsed_time > next_time:
            with open(result_file, 'wb') as fh:
                pickle.dump(dResults,fh)
            next_time = next_time + time_interval

        iS = sequence.Sequence(seq = str(iS_pyfaidx), id = iS_pyfaidx.name)

        log.debug("Work on sequence {}".format(iS))
        ### 1. annotate_de_novo()
        denovo_repeat_list = iS.detect(denovo = True, repeat = {"calc_pValue": True})
        log.debug(denovo_repeat_list.repeats)
        for iTR in denovo_repeat_list.repeats:
            iTR.model = None

        ### 2. annotate_TRs_from_hmmer()
        lHMM = dHMM_annotation[iS.id]
        infoNRuns = len(lHMM)
        log.debug("{} Viterbi runs need to be performed.".format(infoNRuns))
        lHMM = set(lHMM)
        infoNHMM = len(lHMM)
        log.debug("These derive from {} independent HMMs.".format(infoNHMM))
         # Load all HMM pickles needed for the particular sequence.
        if infoNRuns >= 1:
            for hmm_ID in lHMM:
                if hmm_ID not in dHMM:
                    dHMM[hmm_ID] = hmm.HMM.create(format = "pickle", file = os.path.join(hmm_dir, hmm_ID + ".pickle"))

            pfam_repeat_list = iS.detect([dHMM[hmm_ID] for hmm_ID in lHMM], repeat = {"calc_pValue": True})
            for iTR, hmm_ID in zip(pfam_repeat_list.repeats, lHMM):
                iTR.model = hmm_ID
                iTR.TRD = "PFAM"

        ### 3. merge_and_basic_filter()
        all_repeat_list = denovo_repeat_list + pfam_repeat_list
        iS.set_repeat_list(all_repeat_list, REPEAT_LIST_TAG)
        iS.set_repeat_list(denovo_repeat_list, DE_NOVO_ALL_TAG)
        iS.set_repeat_list(pfam_repeat_list, PFAM_ALL_TAG)


        rl_tmp = iS.dRepeat_list[REPEAT_LIST_TAG]
        if iS.dRepeat_list[REPEAT_LIST_TAG]:
            for iB in basic_filter.values():
                rl_tmp = rl_tmp.filter(**iB)
        else:
            rl_tmp = iS.dRepeat_list[REPEAT_LIST_TAG]
        iS.set_repeat_list(rl_tmp, basic_filter_tag)
        iS.set_repeat_list(rl_tmp.intersection(iS.dRepeat_list[PFAM_ALL_TAG]), PFAM_TAG)
        iS.set_repeat_list(rl_tmp.intersection(iS.dRepeat_list[DE_NOVO_ALL_TAG]), DE_NOVO_TAG)

        ### 4. calculate_overlap()

        # Perform common ancestry overlap filter and keep PFAMs
        criterion_pfam_fixed = {"func_name": "none_overlapping_fixed_repeats", "rl_fixed": iS.dRepeat_list[PFAM_TAG], "overlap_type": "common_ancestry"}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**criterion_pfam_fixed)

        # Choose only the most convincing de novo TRs
        criterion_filter_order = {"func_name": "none_overlapping", "overlap": ("common_ancestry", None), "lCriterion": [("pValue", "phylo_gap01"), ("divergence", "phylo_gap01")]}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**criterion_filter_order)

        ### 5. refine_denovo()
        denovo_final = []
        denovo_refined = [None] * len(iS.dRepeat_list[DE_NOVO_ALL_TAG].repeats)
        for i,iTR in enumerate(iS.dRepeat_list[DE_NOVO_ALL_TAG].repeats):
            if not iTR in iS.dRepeat_list[DE_NOVO_TAG].repeats:
                continue
            # Create HMM from TR
            denovo_hmm = hmm.HMM.create(format = 'repeat', repeat = iTR)
            # Run HMM on sequence
            denovo_refined_rl = iS.detect(lHMM = [denovo_hmm])
            append_refined = False
            if denovo_refined_rl and denovo_refined_rl.repeats:
                iTR_refined = denovo_refined_rl.repeats[0]
                iTR_refined.TRD = iTR.TRD
                iTR_refined.model = "cpHMM"
                denovo_refined[i] = iTR_refined
                # Check whether new and old TR overlap. Check whether new TR is significant. If not both, put unrefined TR into final.
                if repeat_list.two_repeats_overlap("shared_char", iTR, iTR_refined):
                    rl_tmp = repeat_list.Repeat_list([iTR_refined])
                    log.debug(iTR_refined.msa)
                    for iB in basic_filter.values():
                        rl_tmp = rl_tmp.filter(**iB)
                    if rl_tmp.repeats:
                        append_refined = True
            else:
                denovo_refined[i] = False
            if append_refined:
                denovo_final.append(iTR_refined)
            else:
                denovo_final.append(iTR)

        iS.set_repeat_list(repeat_list.Repeat_list(denovo_refined), DE_NOVO_REFINED_TAG)
        iS.set_repeat_list(repeat_list.Repeat_list(denovo_final), DE_NOVO_FINAL_TAG)
        iS.set_repeat_list(iS.dRepeat_list[DE_NOVO_FINAL_TAG] + iS.dRepeat_list[PFAM_TAG], FINAL_TAG)

        dResults[iS.id] = iS

    ### 6.a Save results as pickle
    with open(result_file, 'wb') as fh:
        pickle.dump(dResults,fh)

    ### 6.b Save serialized results
    with open(result_file_serialized, 'w') as fh_o:

        if format == 'tsv':
            header = ["ID", "MSA", "begin", "pValue", "lD", "n", "nD", "TRD", "model"]
        fh_o.write("\t".join(header))

        for iS in dResults.values():
            for iTR in iS.dRepeat_list[FINAL_TAG].repeats:
                if format == 'tsv':
                    try:
                        data = [str(i) for i in [iS.id, " ".join(iTR.msa), iTR.begin, iTR.pValue("phylo_gap01"), iTR.lD, iTR.n, iTR.nD, iTR.TRD, iTR.model]]
                    except:
                        print(iTR)
                        raise Exception("(Could not save data for the above TR.)")

                fh_o.write("\n" + "\t".join(data))




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
                iTR.TRD = "PFAM"
            dTR[iS.id] = iRL
        else:
           dTR[iS.id] = None

    with open(result_file, 'wb') as fh:
        pickle.dump(dTR, fh)

    print("\n" + "\t".join(["infoNRuns", "infoNHMM"]))
    print("\n" + "\t".join([str(infoNRuns), str(infoNHMM)]))
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
            iTR.model = None
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

    basic_filter = config['filter']['basic']['dict']
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
            for iB in basic_filter.values():
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
        criterion_pfam_fixed = {"func_name": "none_overlapping_fixed_repeats", "rl_fixed": iS.dRepeat_list[PFAM_TAG], "overlap_type": "common_ancestry"}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**criterion_pfam_fixed)

        # Choose only the most convincing de novo TRs
        criterion_filter_order = {"func_name": "none_overlapping", "overlap": ("common_ancestry", None), "lCriterion": [("pValue", "phylo_gap01"), ("divergence", "phylo_gap01")]}
        iS.dRepeat_list[DE_NOVO_TAG] = iS.dRepeat_list[DE_NOVO_TAG].filter(**criterion_filter_order)

    with open(result_file, 'wb') as fh:
        pickle.dump(lSequence, fh)

    print("DONE")


def refine_denovo(sequences_file, result_file):
    ''' Refine denovo TRs.

    Refine denovo TRs from the DE_NOVO_TAG ``repeat_list``.
    If the refined de novo TR
        - exists
    append it to the DE_NOVO_REFINED_TAG ``repeat_list``. Otherwise, append False. Append
    None for TRs that are in DE_NOVO_ALL_TAG, but not in DE_NOVO_TAG, such that finally,
    all TRs in DE_NOVO_ALL_TAG have a corresponding entry (which might be None or False)
    in DE_NOVO_REFINED_TAG.

    If the refined de novo TR
        - exists
        - overlaps with the original de novo TR
        - passes the basic filtering test
    append it to the FINAL_TAG ``repeat_list``. Otherwise, append the original de novo TR.

     Args:
         sequences_file (str): Path to the pickle file containing a list of ``Sequence``
            instances.
         result_file (str): Path to the result file.

     Raises:
        Exception: If the pickle ``sequences_file`` cannot be loaded
    '''

    basic_filter = config['filter']['basic']['dict']

    try:
        with open(sequences_file, 'rb') as fh:
            lSequence = pickle.load(fh)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    for iS in lSequence:
        log.debug(iS.id)
        denovo_final = []
        denovo_refined = [None] * len(iS.dRepeat_list[DE_NOVO_ALL_TAG].repeats)
        for i,iTR in enumerate(iS.dRepeat_list[DE_NOVO_ALL_TAG].repeats):
            if not iTR in iS.dRepeat_list[DE_NOVO_TAG].repeats:
                continue
            # Create HMM from TR
            denovo_hmm = hmm.HMM.create(format = 'repeat', repeat = iTR)
            # Run HMM on sequence
            denovo_refined_rl = iS.detect(lHMM = [denovo_hmm])
            append_refined = False
            if denovo_refined_rl and denovo_refined_rl.repeats:
                iTR_refined = denovo_refined_rl.repeats[0]
                iTR_refined.TRD = iTR.TRD
                iTR_refined.model = "cpHMM"
                denovo_refined[i] = iTR_refined
                # Check whether new and old TR overlap. Check whether new TR is significant. If not both, put unrefined TR into final.
                if repeat_list.two_repeats_overlap("shared_char", iTR, iTR_refined):
                    rl_tmp = repeat_list.Repeat_list([iTR_refined])
                    log.debug(iTR_refined.msa)
                    for iB in basic_filter.values():
                        rl_tmp = rl_tmp.filter(**iB)
                    if rl_tmp.repeats:
                        append_refined = True
            else:
                denovo_refined[i] = False
            if append_refined:
                denovo_final.append(iTR_refined)
            else:
                denovo_final.append(iTR)

        iS.set_repeat_list(repeat_list.Repeat_list(denovo_refined), DE_NOVO_REFINED_TAG)
        iS.set_repeat_list(repeat_list.Repeat_list(denovo_final), DE_NOVO_FINAL_TAG)
        iS.set_repeat_list(iS.dRepeat_list[DE_NOVO_FINAL_TAG] + iS.dRepeat_list[PFAM_TAG], FINAL_TAG)

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
        Exception: If ``format`` is not implemented.
    '''

    if not format in ['tsv']:
        raise Exception("Serializiation format not implemented: {}".format(format))

    lFiles = [file for file in os.listdir(sequences_dir) if file.endswith(".pickle")]

    with open(result_file, 'w') as fh_o:
        for iFile in lFiles:
            sequences_file = os.path.join(sequences_dir, iFile)
            try:
                with open(sequences_file, 'rb') as fh_i:
                    lSequence = pickle.load(fh_i)
            except:
                raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

            if format == 'tsv':
                header = ["ID", "MSA", "begin", "pValue", "lD", "n", "nD", "TRD", "model"]
            fh_o.write("\t".join(header))

            for iS in lSequence:
                for iTR in iS.dRepeat_list[FINAL_TAG].repeats:
                    if format == 'tsv':
                        try:
                            data = [str(i) for i in [iS.id, " ".join(iTR.msa), iTR.begin, iTR.pValue("phylo_gap01"), iTR.lD, iTR.n, iTR.nD, iTR.TRD, iTR.model]]
                        except:
                            print(iTR)
                            raise Exception("(Could not save data for the above TR.)")

                    fh_o.write("\n" + "\t".join(data))

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
        annotate_de_novo(pars["input"], pars["output"], detector = pars['detector'])
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
        serialize_annotations(pars["input"], pars["output"], pars["format"])

def read_commandline_arguments():
    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('method', metavar='method_name', type=str,
                       help='The name of the method to be executed.')
    parser.add_argument('-i', '--input', type=str, required=True,
                       help='The path to the input file.')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='The path to the output file.')
    parser.add_argument('-seq','--sequence_type', type=str,
                       help='The sequence type: -seq AA or -seq DNA')
    parser.add_argument('-rep','--repeat_files', nargs='+', type=str,
                       help='The repeat files')
    parser.add_argument('-f','--format', type=str,
                       help='The serialization format, e.g. tsv')
    parser.add_argument('-ov','--overlap_type', nargs='+', type=str,
                       help='The overlap type, e.g. "common_ancestry shared_char"')
    parser.add_argument('-d', '--detector', type=str,
                        help='The de novo tandem repeat detectors. For example: -d T-REKS')
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