#!/usr/bin/python
"""
Filter search_hmm results according to various statistics

May be run as a module, e.g. `python -m tral.search.filter_hmm -h`

Filtering is performed by `filter_search_results`. Other methods provide
I/O and supporting roles.

@author Spencer Bliven <sbliven@ucsd.edu>
"""

import argparse
import logging
import gzip
import random
import itertools
import re

from ..hmm import hmm, hmm_io, hmm_viterbi
from ..sequence import sequence
from Bio import SeqIO
from .search_hmm import TralHit, opengzip


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


def filter_search_results(results, repeats=2, log_odds=8.0):
    """Get the set of hits passing the specified filter

    :param results: TSV file with TRAL hits
    :param repeats: minimum number of repeats
    :param log_odds: minimum log odds score
    :return: Generator for all hits
    """
    hits = parse_hits(results)
    filtered = filter(lambda hit: hit.logodds >= log_odds and \
                                  hit.states and \
                                  count_repeats(hit.states) >= repeats,
                      hits)
    for hit in filtered:
        yield hit


def filter_fasta(databasefile, outfile, hits, usedescription=True):
    """Filter a fasta file to IDs contained in the hits

    :param databasefile: Name of input fasta. May be gzipped
    :param outfile: Name of output fasta. Uncompressed.
    :param hits: Iterable of TralHits
    :param usedescription: Include full description (header line) from the input
        fasta database. If false, only the name (first word after the '>') will
        be used.
    """
    ids = {hit.id for hit in hits}
    with opengzip(databasefile) as database:
        with open(outfile, 'w') as results:
            logging.info("Loading database from %s", databasefile)
            seqs = SeqIO.parse(database, 'fasta')
            filtered = filter(lambda seq: seq.name in ids, seqs)
            if not usedescription:
                def replacedescription(seq):
                    seq.description = seq.name
                    return seq
                filtered = map(replacedescription, filtered)
            SeqIO.write(filtered, results, 'fasta')



def parse_hits(results):
    """Generator for hits from a search_hmm output TSV file
    """
    # Filename
    if not hasattr(results, 'readlines'):
        with open(results, 'r') as f:
            for hit in parse_hits(f):
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

def write_hits(hits, outfile):
    """Write a collection of hits to a TSV file

    :param hits: iterable of TralHit
    :param outfile: filename
    """
    with open(outfile, 'w') as results:
        # header
        results.write(TralHit.header)
        results.write("\n")

        for hit in hits:
            results.write(hit.to_line())
            results.write("\n")

def write_treks(databasefile, outfile, hits, hmm=None):
    """Write a collection of hits as a TREKS file

    :param hits: iterable of TralHit
    :param outfile: filename
    """
    with opengzip(databasefile) as database:
        with open(outfile, 'w') as results:
            logging.info("Loading database from %s", databasefile)
            seqs = {rec.name: rec for rec in SeqIO.parse(database, 'fasta')}

            for hit in hits:
                rec = seqs[hit.id]
                results.write(hit.to_treks(rec.seq, hmm=hmm))
                results.write("\n")

def match_seqs(hits, fastafile):
    """Pair a series of hits with sequences

    Requires loading fastafile into memory

    Return: generator of (TralHit, SeqRecord) tuples
    """
    with opengzip(fastafile) as database:
        seqs = {rec.name: rec for rec in SeqIO.parse(database, 'fasta')}
        for hit in hits:
            seq = seqs[hit.id]
            yield (hit, seq)

def main(args=None):
    "filter_hmm main method"
    parser = argparse.ArgumentParser(description='Filter hits from search_hmm')
    parser.add_argument("hits", help="TSV file, as produced by search_hmm, containing HMM hits")
    parser.add_argument("database", help="database to filter, in fasta format (may be gzip compressed)")
    group = parser.add_argument_group('outputs')
    group.add_argument("-f", "--filtered-fasta", help="output fasta file, filtered by hits")
    group.add_argument("-o", "--filtered-tsv", help="Filtered TSV file")
    group.add_argument("-x", "--filtered-treks", help="Filtered T-Reks file")
    group = parser.add_argument_group('Thresholds')
    group.add_argument("-r", "--min-repeats", help="Minimum number of repeats", type=float, default=2.0)
    group.add_argument("-t", "--log-odds", help="Threshold for minimum log-odds ratio", type=float, default=8.0)
    parser.add_argument("--preserve-header", help="Include the full header in FASTA output. "
        "Otherwise, just the identifier is used to match TSV and TREKS.",
        default=False, action="store_true")
    parser.add_argument("-v", "--verbose", help="Long messages",
                        default=False, action="store_true")
    args = parser.parse_args(args)

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO if args.verbose else logging.WARN)

    # parse and filter hits
    hits = filter_search_results(args.hits, args.min_repeats, args.log_odds)

    # preserve iterable if multiple output formats
    if bool(args.filtered_tsv) + bool(args.filtered_fasta) + bool(args.filtered_treks) > 1:
        hits = list(hits)

    # output filtered tsv
    if args.filtered_tsv:
        write_hits(hits, args.filtered_tsv)
    # output filtered fasta
    if args.filtered_fasta:
        filter_fasta(args.database, args.filtered_fasta, hits,
                     usedescription=args.preserve_header)
    # output filtered treks
    if args.filtered_treks:
        write_treks(args.database, args.filtered_treks, hits, hmm=None)

if __name__ == "__main__":
    main()