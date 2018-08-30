#!/usr/bin/python
"""
@author Spencer Bliven <sbliven@ucsd.edu>
"""

import argparse
import logging
import gzip
import random
import itertools
from io import StringIO
import re

from tral.hmm import hmm, hmm_io, hmm_viterbi
from tral.sequence import sequence
from Bio import SeqIO


class TralHit(object):
    """Encapsulates a search result"""
    def __init__(self, id, prob, logodds, states):
        self.id=id
        self.prob=prob
        self.logodds=logodds
        self.states=states

    header = "ID\tProb\tLO\tStates"

    @classmethod
    def parse_line(Cls, line):
        fields = line.strip("\n").split('\t')

        if len(fields) != 4:
            raise ValueError("Illegal format")

        id = fields[0]
        prob = float(fields[1])
        lo = float(fields[2])
        states = fields[3].split(' ') if fields[3] else []
        return Cls(id, prob, lo, states)

    def to_line(self):
        """Converts this hit to a tab-delimited line

        Does not include endline.

        May be recreated by `parse_line(line)`

        Return: (str)
        """
        return "\t".join([
                self.id,
                str(self.prob),
                str(self.logodds),
                " ".join(self.states) if self.states else ""])

    def __repr__(self):
        return f"{self.__class__.__name__}({self.id!r},{self.prob!r}, {self.logodds!r}, {len(self.states)} states)"

    def _to_align_matrix(self, statelist, seq, ignoredstates=["N", "C"]):
        """Constructs a 2D matrix of the repeat alignment from the states of this hit.

        Args:
            - statelist (list of string): a sorted list of all states in this alignment
            - seq (string): sequence for this hit. Should match self.states in length
            - ignoredstates (list of string): list of non-aligned states

        Returns: (3D list of strings)
            The result is a rectangular matrix, where the first index is the repeat number,
            the second index is the column (with the same length as statelist), and the element
            is a list of zero or more letters from seq.
        """
        # reverse index
        statepos = {state: i for i, state in enumerate(statelist)}

        assert len(seq) == len(self.states)

        align = []

        row = []
        for state, res in zip(self.states, seq):
            if state in ignoredstates:
                continue
            assert state in statelist  # should be the only other states

            pos = statepos[state]
            if pos < len(row) - 1:
                # new row
                row.extend([] for i in range(len(statelist) - len(row)))
                align.append(row)
                row = []

            row.extend([] for i in range(pos + 1 - len(row)))

            row[-1].append(res)
        row.extend([] for i in range(len(statelist) - len(row)))
        align.append(row)
        return align

    def _to_repeat_alignment(self, align, statelist, stateheader=False):
        """Constructs a multiple alignment string for this hit.

        Gaps are indicated by -;

        Args:
            - align (3D list of strings): output of self._to_align_matrix
            - statesheader (bool): Add header rows giving the state names for each column
        Returns: (3D list of strings)

        """
        # Find max width of each column
        widths = [max(len(align[row][col]) for row in range(len(align))) for col in range(len(statelist))]

        # Convert to strings
        alignstr = StringIO()

        if stateheader:
            maxlen = max(len(s) for s in statelist)
            for row in range(maxlen):
                for col in range(len(statelist)):
                    w = widths[col]
                    char = statelist[col][row:(row+1)]
                    alignstr.write((char if char else " ")*w)
                alignstr.write("\n")
        for row in range(len(align)):
            for col in range(len(statelist)):
                w = widths[col]
                alignstr.write("{:-<{}}".format("".join(align[row][col]), w))
            alignstr.write("\n")

        return alignstr.getvalue()

    def to_treks(self, sequence, hmm=None):
        """Converts this hit to T-REKS format

        Args:
            sequence (str or Bio.SeqRecord or Bio.Seq): sequence for this hit
            hmm (HMM): (optional) the HMM used to generate this hit. If specified,
                results in more accurate set of states
        Return: (str)

        See: tral.sequence.repeat_detection_io.treks_get_repeats
        """
        # Example:
        # >sp|Q10567|AP1B1_HUMAN AP-1 complex subunit beta-1 OS=Homo sapiens OX=9606 GN=AP1B1 PE=1 SV=2
        # Length: 4 residues - nb: 5  from  625 to 645 - Psim:0.7 region Length:21
        # G-D-LL
        # G-D-LL
        # NLD-L-
        # G-PPVS
        # G-P-PL
        # **********************


        # Handle case without any repeats
        if not self.states:
            return "repeat not found in sequence {}".format(self.id)

        if hasattr(sequence,"seq"):
            seq = sequence.seq
        else:
            seq = sequence
        if len(seq) != len(self.states):
            raise ValueError("Mismatched sequence lengths for " + self.id)

        out = StringIO()

        # identifier
        desc = sequence.description if hasattr(sequence,"description") and sequence.description else self.id
        out.write(">")
        out.write(desc)
        out.write("\n")

        # Compute header info
        usedstates = set(self.states)

        # Guess HMM length
        if hmm:
            hmm_len = hmm.l_effective
            states = [state for m_i in zip(hmm.match_states, hmm.insertion_states)
                      for state in m_i
                      if state in usedstates]
            ignored = hmm.terminal_states
        else:
            # Check that state names are expected
            statepattern = re.compile(r"[NC]|[MI]\d+")
            if not all(statepattern.match(s) for s in usedstates):
                raise Exception("Unknown state naming pattern. HMM required.")

            # Guess based on state names
            hmm_len = max(int(a[1:]) for a in self.states if len(a) > 1)
            states = [state for i in range(1, hmm_len+1)
                      for state in ("M%d" % i, "I%d" % i)
                      if state in usedstates]
            ignored = ["N", "C"]

        # Compute alignment
        align = self._to_align_matrix(states, seq, ignored)

        # Compute bounds. Should be 0-based inclusive
        start = min(i for i, s in enumerate(self.states) if s not in ignored)
        end = max(i for i, s in enumerate(self.states) if s not in ignored)

        out.write("Length: {} residues - ".format(hmm_len))
        out.write("nb: {}  ".format(len(align)))
        out.write("from  {} to {} - ".format(start, end))
        out.write("Psim:{:.2f} ".format(self.prob))
        out.write("region Length:{}".format(end - start + 1))
        out.write("\n")

        # write alignment
        out.write(self._to_repeat_alignment(align, states))
        out.write("**********************\n")

        return out.getvalue()


def opengzip(filename):
    """Open a file, which may optionally be gzip'd"""
    magic = None
    with open(filename, 'br') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':  # gzip
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'rt')


def shuffle_seq(seq):
    """Randomly shuffle a sequence

    Args:
        - seq (Bio.Seq): input sequence

    Returns:
        Bio.Seq.MutableSeq with suffled characters

    """
    mutable = seq.tomutable()
    random.shuffle(mutable)
    return mutable


def tralsearch(hmmfile, databasefile, outfile, start=0, n=0, shuffle=False):
    circular_profile = hmm.HMM.create(input_format='hmmer', file=hmmfile)

    with opengzip(databasefile) as database:
        with open(outfile, 'w') as results:
            logging.info("Loading database from %s", databasefile)
            # Use BioPython over TRAL for performance
            # seqs = sequence.Sequence.create(database, input_format='fasta')
            # logging.info("Loaded %d sequences", len(list(seqs)))
            seqs = SeqIO.parse(database, 'fasta')

            # header
            results.write(TralHit.header)
            results.write("\n")

            for record in itertools.islice(seqs, start, start + n if n > 0 else None):

                seq = record.seq
                if shuffle:
                    seq = shuffle_seq(seq)

                hmm_results, prob = hmm_viterbi.viterbi_with_prob(circular_profile, seq)

                odds = hmm_viterbi.logodds(circular_profile, seq, prob)

                hit = TralHit(record.name, prob, odds, hmm_results)
                results.write(hit.to_line())
                results.write("\n")
                results.flush()


def main(args=None):
    parser = argparse.ArgumentParser(description='Align a cpHMM against a database using TRAL')
    parser.add_argument("hmm", help="cpHMM filename")
    parser.add_argument("database", help="database to search, in fasta format (may be gzip compressed)")
    parser.add_argument("results", help="results file")
    parser.add_argument("-v", "--verbose", help="Long messages",
                        default=False, action="store_true")
    parser.add_argument("-s", "--start", help="Index of database sequence to start with (default: 0)",
                        type=int, default=0)
    parser.add_argument("-n", "--dbsize", help="Number of database sequences to include (default: 0 to include all)",
                        type=int, default=0)
    parser.add_argument("-r", "--shuffle", help="Shuffle database sequences",
                        default=False, action="store_true")
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO if args.verbose else logging.WARN)

    tralsearch(args.hmm, args.database, args.results, start=args.start, n=args.dbsize, shuffle=args.shuffle)

if __name__ == "__main__":
    main()