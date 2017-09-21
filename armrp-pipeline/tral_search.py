#!/usr/bin/python
"""
@author Spencer Bliven <sbliven@ucsd.edu>
"""

import sys
import os
import argparse
import logging
import gzip
import itertools
from tral.hmm import hmm, hmm_io, hmm_viterbi
from tral.sequence import sequence
from Bio import SeqIO


def opengzip(filename):
    """Open a file, which may optionally be gzip'd"""
    magic = None
    with open(filename, 'br') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':  # gzip
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'rt')


def tralsearch(hmmfile, databasefile, outfile, start=0, n=0):
    circular_profile = hmm.HMM.create(input_format='hmmer', file=hmmfile)

    with opengzip(databasefile) as database:
        with open(outfile, 'w') as results:
            logging.info("Loading database from %s", databasefile)
            # Use BioPython over TRAL for performance
            # seqs = sequence.Sequence.create(database, input_format='fasta')
            # logging.info("Loaded %d sequences", len(list(seqs)))
            seqs = SeqIO.parse(database, 'fasta')

            for seq in itertools.islice(seqs, start, start + n if n > 0 else None):
                hmm_results, prob = hmm_viterbi.viterbi_with_prob(circular_profile, seq.seq)

                results.write(seq.name)
                results.write("\t")
                results.write(str(prob))
                results.write("\t")
                results.write(" ".join(hmm_results))
                results.write("\n")
                results.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Align a cpHMM against a database using TRAL')
    parser.add_argument("hmm", help="cpHMM filename")
    parser.add_argument("database", help="database to search, in fasta format")
    parser.add_argument("results", help="results file")
    parser.add_argument("-v", "--verbose", help="Long messages",
        dest="verbose", default=False, action="store_true")
    parser.add_argument("-s", "--start", help="Index of database sequence to start with (default: 0)",
                        type=int, default=0)
    parser.add_argument("-n", "--dbsize", help="Number of database sequences to include (default: 0 to include all)",
                        type=int, default=0)
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO if args.verbose else logging.WARN)

    tralsearch(args.hmm, args.database, args.results, start=args.start, n=args.dbsize)
