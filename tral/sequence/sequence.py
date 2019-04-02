# (C) 2015 Elke Schaper

"""
    :synopsis: The Sequence Class.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import pickle
import re

from tral import configuration
from tral.repeat import repeat, repeat_align
from tral.repeat_list import repeat_list
from tral.hmm import hmm, hmm_viterbi
from tral.sequence import repeat_detection_run, sequence_io

CONFIG = configuration.Configuration.instance().config
LOG = logging.getLogger(__name__)


class Sequence:

    """ A ``Sequence`` describes either a protein or a DNA sequence.

    ``Sequence`` contains methods that act on single sequences, for example:

    *   tandem repeats

    Attributes:
        seq (str): The sequence.
        seq_standard_aa (str): The sequence with standard amino acids only
    """

    def __init__(self, seq, name=None, sequence_type="AA"):

        if not isinstance(seq, str):
            raise Exception('The seq value is not a String')
        self.seq = seq.upper()
        self.sequence_type = sequence_type
        for i in self.seq:
            if i not in CONFIG[self.sequence_type]['all_chars']:
                raise Exception("{} is not in CONFIG[{}]['all_chars']: {}"
                                .format(i, self.sequence_type, CONFIG['all_chars']))

        self.seq_standard_aa = repeat.standardize(self.seq, self.sequence_type)

        if name:
            self.name = name

        self.d_annotations = {}
        self.d_repeatlist = {}

    def create(file, input_format):
        """ Create sequence(s) from file.

        Create sequence(s) from file.

        Args:
            file (str): Path to input file
            format (str):  Either "fasta" or "pickle"

        .. todo:: Write checks for ``format`` and ``file``.
        """

        if input_format == 'fasta':
            l_seq = sequence_io.read_fasta(file)
            return [Sequence(iSeq, iID) for iSeq, iID in l_seq]
        if input_format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        else:
            raise Exception("Input format {} is not implemented for"
                            "sequence.create()".format(input_format))

    def write(self, file, file_format):
        """ Write sequence to file.

        Write sequence to file using one of two formats.

        Args:
            file (str): Path to output file
            format (str):  Either "fasta" or "pickle"

        .. todo:: Write checks for ``format`` and ``file``.

        """

        if file_format == 'fasta':
            sequence_io.write(self.seq, file)
        elif file_format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        else:
            raise Exception("Output format {} is not implemented for",
                            "sequence.write()".format(file_format))

    def detect(self, lHMM=None, denovo=None, aligner='mafft', sequence_type='AA', **kwargs):
        """ Detects tandem repeats on ``self.seq`` from 2 possible sources.

        A list of ``Repeat`` instances is created for tandem repeat detections
        on the sequence from two possible sources:

        * Sequence profile hidden Markov models ``HMM``
        * de novo detection algorithms.

        Args:
            hmm (HMM): A list of ``HMM`` instances.
            denovo (bool): boolean
            *kwargs: Parameters fed to denovo TR prediction and/or Repeat
                instantiation. E.g. ``repeat = {"calc_score": True}``

        Returns:
            A ``RepeatList`` instance
        """

        if lHMM:
            if not isinstance(lHMM, list):
                raise Exception('The lHMM value is not a list.')
            for iHMM in lHMM:
                if not isinstance(iHMM, hmm.HMM):
                    raise Exception('At least one list element in the lHMM'
                                    'value is not a valid instance of the HMM'
                                    'class.')

            repeats = []
            for iHMM in lHMM:
                # Detect TRs on self.seq with hmm using the Viterbi algorithm.
                most_likely_path = iHMM.viterbi(self.seq)
                LOG.debug("most_likely_path: {}".format(most_likely_path))
                if not most_likely_path:
                    continue
                unaligned_msa = hmm_viterbi.hmm_path_to_non_aligned_tandem_repeat_units(
                    self.seq,
                    most_likely_path,
                    iHMM.l_effective)
                if len(unaligned_msa) > 1:
                    # Align the msa
                    aligned_msa = repeat_align.realign_repeat(unaligned_msa, aligner, sequence_type)
                    if len(aligned_msa) > 1:
                        # Create a Repeat() class with the new msa
                        if 'repeat' in kwargs:
                            repeats.append(repeat.Repeat(aligned_msa,
                                                         **kwargs['repeat']))
                        else:
                            repeats.append(repeat.Repeat(aligned_msa))

            # Set begin coordinate for all repeats
            for i_repeat in repeats:
                self.repeat_in_sequence(i_repeat)

            return repeat_list.RepeatList(repeats)

        elif lHMM == []:
            LOG.debug("lHMM == []")
            return None

        elif denovo:
            if 'detection' in kwargs:
                predicted_repeats = repeat_detection_run.run_detector(
                    [self],
                    **kwargs['detection'])[0]
            else:
                predicted_repeats = repeat_detection_run.run_detector([self])[0]

            LOG.debug("predicted_repeats: {}".format(predicted_repeats))
            repeats = []

            for jTRD, jlTR in predicted_repeats.items():
                for iTR in jlTR:
                    if len(iTR.msa) < 3 or len(iTR.msa[0]) < 10:
                        #print("Tandem repeats with less than 3 units or unit length of less than 10 characters will not be realigned.")
                        pass
                    else:
                        # realign tandem repeats   
                        # print("Starting of realignment.")                     
                        iTR.msa = repeat_align.realign_repeat(iTR.msa, aligner, sequence_type)
                        # TODO: include option to set gamma distribution

                    if 'repeat' in kwargs:
                        iTR = repeat.Repeat(iTR.msa, begin=iTR.begin,
                                            **kwargs['repeat'])
                    else:
                        iTR = repeat.Repeat(iTR.msa, begin=iTR.begin)

                    # Consider only tandem repeats that have a repeat unit
                    # predicted to be at least one character long.
                    if iTR.l_effective > 0:

                        # Save l, n, MSA, TRD, scores, sequence_type, position
                        # in sequence of given type
                        iTR.TRD = jTRD

                        # Sanity check repeat and set begin coordinate for
                        # all repeats
                        if not self.repeat_in_sequence(iTR):
                            LOG.debug("The tandem repeat is not part of"
                                      "the sequence. Detector: %s", iTR.TRD)
                            continue

                        repeats.append(iTR)

            return repeat_list.RepeatList(repeats)

        else:
            raise Exception("Either require denovo detection, or provide an",
                            "HMM")

    def get_repeatlist(self, tag):
        """ Retrieve `repeatlist` from this `sequence` instance.

        Retrieve `repeatlist` from this `sequence` instance. Access
        `repeatlist` as self.d_repeatlist[tag]

        Args:
            tag (str): A identifier for the repeat_list

        Returns:
            RepeatList: A repeat_list instance.
        """

        try:
            return self.d_repeatlist[tag]
        except ValueError:
            logging.error("RepeatList %s not in Sequence", tag)

    def set_repeatlist(self, repeatlist, tag):
        """ Add `repeatlist` as attribute to this `sequence` instance.

        Add `repeatlist` as attribute to this `sequence` instance. Access
        `repeatlist` as self.d_repeatlist[tag]

        Args:
            repeatlist (RepeatList): A repeat_list instance.
            tag (str): A identifier for the repeat_list
        """

        self.d_repeatlist[tag] = repeatlist

    def annotate(self, data, tag):

        self.d_annotations[tag] = data

    def get_annotation(self, tag):

        if tag in self.d_annotations:
            return self.d_annotations[tag]
        else:
            return []

    def repeat_in_sequence(self, my_repeat):
        """ Sanity check whether the `repeat` is part of this `sequence`. In
        case, calculate the position of the `repeat` within the `sequence`.

        If yes: Return True, set repeat.begin to corrected value if necessary.
        If no: Return False.
        Perform sanity check on sequences where all amino acids are or are
        converted to standard amino acids.

        Args:
            sequence (sequence): A sequence instance.

        Returns:
            bool: True if repeat is part of sequence, else false

        .. todo:: Decide whether save_original_msa is needed here.
        """

        repeat_sequence = repeat.get_repeat_sequence(my_repeat.msa_standard_aa)
        # The first letter in the sequence is counted as 1
        # (not 0, as in Python):
        starts = [m.start() + 1 for m in re.finditer(repeat_sequence,
                                                     self.seq_standard_aa)]

        if len(starts) != 0:  # Is the tandem repeat predicted correctly?
            if not hasattr(my_repeat, "begin") or my_repeat.begin not in starts:
                my_repeat.begin = starts[0]
            my_repeat.save_original_msa(self.seq)
            return True
        else:
            return False
