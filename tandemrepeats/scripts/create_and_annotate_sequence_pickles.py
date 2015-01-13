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


def create_and_annotate_seq_pickles(sequence_dir, output_path, annotation_file = None, lFiles = None):

    ''' Create ``Sequence`` instances from fasta files and annotate with data.

    Create ``Sequence`` instances from fasta files and annotate with data.

    Args:
         sequence_dir (str): Path to the dir containing *.fasta sequence files.
         output_path (str): Path to directory where the output files are saved.
         annotation_file (str): Path to a tab separated file with annotations to the
            sequences.

    Raises:
        Exception: If the pickle ``annotation_file`` cannot be loaded.
    '''

    if annotation_file:
        try:
            with open(annotation_file, 'rb') as fh:
                annotations = pickle.load(fh)
        except:
            raise Exception("Cannot load sequence annotation file annotation_file: {}".format(annotation_file))

    if not lFiles:
        lFiles = list(os.listdir(sequence_dir))
        lFiles = [file for file in lFiles if file.endswith(".fasta")]

    for file in lFiles:
        lSeq = sequence.Sequence.create(file = os.path.join(sequence_dir, file), format = 'fasta')

        if annotation_data_file:
            for iS in lSeq:
                # The fasta files sequence ID follow the pattern "sp|SPID|SPNAME"
                # Here, we extract SPID
                iS.id_long = iS.id
                iS.id = iS.id.split("|")[1]
                if iS.id in annotations.keys():
                    print(annotations[iS.id])
                    iS.annotate(annotations[iS.id], "PFAM")

        output_file = os.path.join(output_path, file.replace("fasta", "pickle"))
        with open(output_file, 'wb') as fh:
            pickle.dump(lSeq, fh)

    print("DONE")


def main():

    pars = read_commandline_arguments()

    kwargs = {"sequence_dir": pars["input"], "output_path": pars["output_path"]}

    if "annotation_file" in pars:
        kwargs["annotation_file"] = pars["annotation_file"]
    if "file_list" in pars:
        kwargs["lFiles"] = pars["file_list"]

    create_and_annotate_seq_pickles(**kwargs)


def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process create hmm pickles options')
    parser.add_argument('-i','--input', type=str, required=True,
                       help='The path to the sequence directory containing .fasta files, e.g. /path/to/sequence_dir')
    parser.add_argument('-o','--output_path', type=str, required=True,
                       help='The path to the output files, e.g. /path/to/output')
    parser.add_argument('-a', '--annotation_file', type=str,
                       help='The path to the annotation data tab-separated file.')
    parser.add_argument('-f', '--file_list', type=str, nargs='+',
                       help='A list of files to work on.')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()
