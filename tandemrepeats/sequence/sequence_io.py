import Bio.Seq
import logging
import os
import re


logger = logging.getLogger('root')

################################## READ SEQUENCE #########################################


def read_fasta(file, indices = None):

    """ Read all sequences from a fasta file.

    Read all sequences from a fasta file.
    At current, the Biopython SeqIO parser is used.

    Args:
        file (str): Path to input file
        start ([int, int]): Index of the first returned sequence, and the first not returned sequence.

    .. todo:: Write checks for ``format`` and ``file``.
    """

    from Bio import SeqIO

    # Making a list out if the generator object might be overhead for huge fastafiles
    count = 0
    for seq_record in SeqIO.parse(file, "fasta"):
        if indices:
            count += 1
            if count < indices[0]:
                continue
            elif count >= indices[1]:
                break
        yield str(seq_record.seq)


############################## WRITE SEQUENCE #############################################

def write(sequence, sequence_file, sequence_id = "sequence_id_not_defined"):
    ''' save the sequence in fasta format in specified <sequence_file> '''

    with open(sequence_file, 'a') as fastafile:
        fastafile.write(">"+str(sequence_id) + '\n')
        fastafile.write(str(sequence) + '\n')
