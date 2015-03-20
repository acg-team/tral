#!python

import argparse
import csv
import datetime
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
from tral.hmm import hmm_io

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


def workflow(sequences_file, hmm_annotation_file, hmm_dir, result_file, result_file_serialized, format, max_time, time_interval = 3600,  next_time = 3600,  **kwargs):

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

    '''

    start = datetime.datetime.now()
    max_time, time_interval, next_time = int(max_time), int(time_interval), int(next_time)

    try:
        lSequence = Fasta(sequences_file)
    except:
        raise Exception("Cannot load putative pickle file sequences_file: {}".format(sequences_file))

    if not os.path.isdir(hmm_dir):
        try:
            os.makedirs(hmm_dir)
        except:
            raise Exception("hmm_dir does not exists and could not be created: {}".format(hmm_dir))

    try:
        with open(hmm_annotation_file, 'rb') as fh:
            dHMM_annotation = pickle.load(fh)
    except:
        raise Exception("Cannot load hmm_annotation_file: {}".format(hmm_annotation_file))

    basic_filter = config['filter']['basic']['dict']
    basic_filter_tag = config['filter']['basic']['tag']

    # Load previous results
    try:
        if not os.path.isdir(os.path.dirname(result_file)):
            os.makedirs(os.path.dirname(result_file))
    except:
        raise Exception("Could not create path to result_file directory: {}".format(os.path.dirname(result_file)))

    try:
        with open(result_file, 'rb') as fh:
            dResults = pickle.load(fh)
    except:
        log.debug("Could not load previous results file - perhaps non existant: {}".format(result_file))
        dResults = {}

    dHMM = {}
    for iS_pyfaidx in lSequence:

        # If sequence is already included in results: continue.
        if iS_pyfaidx.name in dResults:
            continue

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
        if iS.id in dHMM_annotation:
            lHMM = dHMM_annotation[iS.id]
            infoNRuns = len(lHMM)
            log.debug("{} Viterbi runs need to be performed.".format(infoNRuns))
            lHMM = set(lHMM)
            infoNHMM = len(lHMM)
            log.debug("These derive from {} independent HMMs.".format(infoNHMM))
             # Load all HMM pickles needed for the particular sequence.
            for hmm_ID in lHMM:
                if hmm_ID not in dHMM:
                    dHMM[hmm_ID] = hmm.HMM.create(format = "pickle", file = os.path.join(hmm_dir, hmm_ID + ".pickle"))

            pfam_repeat_list = iS.detect([dHMM[hmm_ID] for hmm_ID in lHMM], repeat = {"calc_pValue": True})
            for iTR, hmm_ID in zip(pfam_repeat_list.repeats, lHMM):
                iTR.model = hmm_ID
                iTR.TRD = "PFAM"
        else:
            pfam_repeat_list = None

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

    print("DONE")


def file_preparation(annotation_data_file, annotation_data_output, hmm_raw_file, hmm_dir):

    create_hmm_files(hmm_raw_file, hmm_dir)

    read_pfam_uniprot(annotation_data_file, annotation_data_output)


def create_hmm_files(hmm_file, output_file):

    for hmmer_probabilities in hmm_io.read(hmm_file):
        iID = hmmer_probabilities['id'].split(".")[0]
        iHMM = hmm.HMM(hmmer_probabilities)
        file = os.path.join(output_file, "{}.pickle".format(iID))
        iHMM.write(file, "pickle")


def read_pfam_uniprot(annotation_data_file, output_file):
    ''' annotation_file from:
        http://www.uniprot.org/uniprot/?query=database:(type:pfam%20AND%20*)&fil=&sort=score '''

    p = {}
    if annotation_data_file:
        try:
            with open(annotation_data_file) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    p[row[0]] = row[1][:-1].split(";")
        except:
            raise Exception("Cannot load sequence annotation file annotation_data_file: {}".format(annotation_data_file))

    with open(output_file, 'wb') as fh:
       pickle.dump(p,fh)


def main():

    pars = read_commandline_arguments()

    # Update configuration
    if "sequence_type" in pars:
        config["sequence_type"] = pars["sequence_type"]
    if "detectors" in pars:
        config["sequence"]["repeat_detection"][config["sequence_type"]] = pars["detectors"]
    if "significance_test" in pars:
        config["repeat_score"]["score_calibration"] = pars["significance_test"]

    if pars["method"] == "file_preparation":
        file_preparation(pars["hmm_annotation_raw"], pars["hmm_annotation"], pars["hmm_raw"], pars["hmm"])
    elif pars["method"] == "workflow":
        workflow(pars["input"], pars["hmm_annotation"], pars["hmm"], pars["output"], pars["output_serialized"], pars["format"], pars["time"])


def read_commandline_arguments():
    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('method', metavar='method_name', type=str,
                       help='The name of the method to be executed.')
    parser.add_argument('-i', '--input', type=str,
                       help='The path to the input file.')
    parser.add_argument('-o', '--output', type=str,
                       help='The path to the output file.')
    parser.add_argument('-os', '--output_serialized',
                       help='The path to the serialized output file.')
    parser.add_argument('-f','--format', type=str,
                       help='The serialization format, e.g. tsv')
    parser.add_argument('-t','--time', type=str,
                       help='The maximum runtime')
    parser.add_argument('-hmm', '--hmm', type=str,
                        help='The path to the dir with HMM pickles')
    parser.add_argument('-hmm_raw', '--hmm_raw', type=str,
                        help='The path to the raw file with HMM models')
    parser.add_argument('-hmm_annotation', '--hmm_annotation', type=str,
                        help='The path the pickle with sequence to HMM mappings')
    parser.add_argument('-hmm_annotation_raw', '--hmm_annotation_raw', type=str,
                        help='The path the raw file with sequence to HMM mappings')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()