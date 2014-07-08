import Bio.Seq
import logging
import os
import re


logger = logging.getLogger('root')

################################## READ SEQUENCE #########################################


def read_from_file(seq_filename, start, number, sequence_type = 'AA'):

    '''Read sequence from fasta. Return a Bio.Seq object '''

    from Bio import SeqIO

    # Making a list out if the generator object might be overhead for huge fastafiles
    for seq_record in list(SeqIO.parse(seq_filename, "fasta"))[start:start+number]:
        seq_record.seq.alphabet = sequence_type
        yield seq_record


############################## WRITE SEQUENCE #############################################

def write(sequence, sequence_file, sequence_id = "sequence_id_not_defined"):
    ''' save the sequence in fasta format in specified <sequence_file> '''

    with open(sequence_file, 'a') as fastafile:
        fastafile.write(">"+str(sequence_id) + '\n')
        fastafile.write(str(sequence) + '\n')
