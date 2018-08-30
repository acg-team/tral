#!/usr/bin/python
"""
@author Spencer Bliven <sbliven@ucsd.edu>
"""

import argparse
import logging
import gzip
import random
import itertools
import re

from tral.hmm import hmm, hmm_io, hmm_viterbi
from tral.sequence import sequence
from Bio import SeqIO
from .tral_search import TralHit, opengzip


def count_repeats(states, hmm_length=0):
    """Count the number of repeats matched by the hmm.

    Partial repeats at the beginning and end of the hmm

    For example, consider a 3-column HMM with the following states:

      N N N M2 M3 M1 I1 M2 M1 M2 C C C

    This would have .6+1+.6=2.3 repeats.

    Args:
        - states (list of str): list of states traversed by the hmm. Only
            match states (starting with 'M') are considered, and are assumed
            to be numbered sequentially from 1
        - hmm_length (int): number of HMM states. If not given, guessed based
            on the observed states (may be inaccurate)


    Returns:
        (float) Number of repeats

    """

    match_re = re.compile('M([0-9]+)', re.IGNORECASE)

    # states are numbered from 1
    first_match = None
    last_match = -1
    circles = 0

    for state in states:
        if state[0] == 'M':  # match state
            # find state index
            match = match_re.search(state)

            if match:
                index = int(match.group(1))
                if not first_match:
                    first_match = index

                # loop around
                if last_match > index:
                    circles += 1

                hmm_length = max(hmm_length, index)

                last_match = index
    if first_match is None:
        return 0

    return -(first_match - 1.) / hmm_length + \
        circles + \
        float(last_match) / hmm_length


def tralfilter(results, repeats=2, log_odds=8.0):
    hits = parsehits(results)
    filtered = filter(lambda hit: hit.logodds >= log_odds and \
                                  hit.states and \
                                  count_repeats(hit.states) >= repeats,
                      hits)
    for hit in filtered:
        yield hit

def filterfasta(databasefile, outfile, hits):
    """Filter a fasta file to IDs contained in the hits

    :param databasefile: Name of input fasta. May be gzipped
    :param outfile: Name of output fasta. Uncompressed.
    :param hits: Iterable of TralHits
    """
    ids = {hit.id for hit in hits}
    with opengzip(databasefile) as database:
        with open(outfile, 'w') as results:
            logging.info("Loading database from %s", databasefile)
            seqs = SeqIO.parse(database, 'fasta')
            filtered = filter(lambda seq: seq.name in ids, seqs)
            SeqIO.write(filtered, results, 'fasta')



def parsehits(results):
    """Generator for hits from a tral_search output TSV file
    """
    # Filename
    if not hasattr(results, 'readlines'):
        with open(results, 'r') as f:
            for hit in parsehits(f):
                yield hit
            return

    # File handle
    lines = iter(results.readlines())
    line = next(lines)
    line = line.strip("\n")
    linenr = 1
    if line is None:
        return  # empty file
    if line != TralHit.header:
        logging.error("Invalid header in results file")
        return

    for line in lines:
        linenr += 1
        try:
            hit = TralHit.parse_line(line)
            yield hit
        except ValueError as e:
            logging.error("{} in results file line {}".format(e, linenr))

def main(args=None):
    parser = argparse.ArgumentParser(description='Filter hits from tral_search')
    parser.add_argument("hits", help="TSV file, as produced by tral_search, containing HMM hits")
    parser.add_argument("database", help="database to filter, in fasta format (may be gzip compressed)")
    parser.add_argument("outfile", help="output fasta file, filtered by hits")
    parser.add_argument("-r", "--min-repeats", help="Minimum number of repeats", type=float, default=2.0)
    parser.add_argument("-t", "--log-odds", help="Threshold for minimum log-odds ratio", type=float, default=8.0)
    parser.add_argument("-v", "--verbose", help="Long messages",
                        default=False, action="store_true")
    args = parser.parse_args(args)

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO if args.verbose else logging.WARN)

    hits = tralfilter(args.hits, args.min_repeats, args.log_odds)
    filterfasta(args.database, args.outfile, hits)

if __name__ == "__main__":
    main()