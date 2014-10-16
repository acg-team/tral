import argparse
from collections import defaultdict
import csv
import logging
import logging.config
import os
import pickle
from Bio import SeqIO

from tandemrepeats.paths import *
from tandemrepeats.sequence import sequence


logging.config.fileConfig(os.path.join(CODEROOT,'tandemrepeats','data','logging.ini'))
log = logging.getLogger('root')


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

    with open('d.pickle', 'wb') as fh:
       pickle.dump(p,fh)


def annotate_seq_pickles(sequence_dir, output_path, annotation_pickle):

    with open(annotation_pickle, 'rb') as fh:
        annotations = pickle.load(fh)

    for file in os.listdir(sequence_dir):
        if file.endswith(".pickle"):
            with open(os.path.join(sequence_dir, file), 'rb') as fh:
                lSeq = pickle.load(fh)

            for iS in lSeq:
                if iS.id in annotations.keys():
                    iS.annotate(annotations[iS.id], "PFAM")

            output_file = os.path.join(output_path, file)
            with open(output_file, 'wb') as fh:
                pickle.dump(lSeq, fh)

    print("DONE")


def create_and_annotate_seq_pickles(sequence_dir, output_path, annotation_data_file = None):

    ''' Create ``Sequence`` instances from fasta files and annotate with data.

    Create ``Sequence`` instances from fasta files and annotate with data.

    Args:
         sequence_dir (str): Path to the dir containing *.fasta sequence files.
         output_path (str): Path to directory where the output files are saved.
         annotation_data_file (str): Path to a tab separated file with annotations to the
            sequences.

    Raises:
        Exception: If the pickle ``annotation_data_file`` cannot be loaded.
    '''

    if annotation_data_file:
        try:
            with open(annotation_data_file) as f:
                reader = csv.reader(f, delimiter="\t")
                d = list(reader)
        except:
            raise Exception("Cannot load sequence annotation file annotation_data_file: {}".format(annotation_data_file))

    for file in os.listdir(sequence_dir):
        if file.endswith(".fasta"):
            lSeq = sequence.Sequence.create(file = os.path.join(sequence_dir, file), format = 'fasta')

            if annotation_data_file:
                for iS in lSeq:
                    lPFAM = [i["accession"] for i in d if d['seq'] == iS['id'] and d['accession_type'] == 'PFAM']
                    iS.annotate(data, "PFAM")

            output_file = os.path.join(output_path, file.replace("fasta", "pickle"))
            with open(output_file, 'wb') as fh:
                pickle.dump(lSeq, fh)

    print("DONE")


def main():

    pars = read_commandline_arguments()
    if "annotation_file" in pars:
        annotate_seq_pickles(pars["sequence_dir"], pars["output_path"], pars["annotation_file"])
    else:
        annotate_seq_pickles(pars["sequence_dir"], pars["output_path"])


def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process create hmm pickles options')
    parser.add_argument('-a', '--annotation_file', type=str,
                       help='The path to the annotation data tab-separated file.')
    parser.add_argument('-s','--sequence_dir', type=str,
                       help='The path to the sequence directory containing .fasta files, e.g. /path/to/sequence_dir')
    parser.add_argument('-o','--output_path', type=str,
                       help='The path to the output files, e.g. /path/to/output')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()
