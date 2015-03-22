# (C) 2015 Elke Schaper

"""
    :synopsis: Input/output for sequences

    .. moduleauthor:: Elke Schaper <elke.schaper@sib-sib.ch>

"""

from Bio import SeqIO
import logging

LOG = logging.getLogger(__name__)

'''  READ SEQUENCE '''


def read_fasta(file, indices=None):
    """ Read all sequences from a fasta file.

    Read all sequences from a fasta file.
    At current, the Biopython SeqIO parser is used.

    Args:
        file (str): Path to input file
        start ([int, int]): Index of the first returned sequence, and the first
                            not returned sequence.

    .. todo:: Write checks for ``format`` and ``file``.
    """

    # Making a list out if the generator object might be overhead for huge
    # fastafiles
    count = 0
    for seq_record in SeqIO.parse(file, "fasta"):
        if indices:
            count += 1
            if count < indices[0]:
                continue
            elif count >= indices[1]:
                break
        yield str(seq_record.seq), seq_record.id


#  WRITE SEQUENCE ##########################################################


def write(sequence, sequence_file, sequence_id="sequence_id_not_defined"):
    """ Write a sequence str to fasta format in specified <sequence_file>

    Write s sequence str to fasta format in specified <sequence_file>

    Args:
        sequence (str): Sequence
        sequence_file (str): Path to the output file
        sequence_id (str): ID of the sequence in the output file.

    """

    with open(sequence_file, 'a') as fastafile:
        fastafile.write(">" + str(sequence_id) + '\n')
        fastafile.write(str(sequence) + '\n')
