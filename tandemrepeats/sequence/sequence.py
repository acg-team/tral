# (C) 2014 Elke Schaper

import logging
import pickle
import os
import re

from tandemrepeats import configuration
from tandemrepeats.repeat import repeat, repeat_align
from tandemrepeats.repeat_list import repeat_list
from tandemrepeats.hmm import hmm, hmm_viterbi
from tandemrepeats.sequence import repeat_detection_run, sequence_io
from tandemrepeats.paths import *

c = configuration.Configuration.Instance()
config = c.config

log = logging.getLogger(__name__)

class Sequence:

    """ A ``Sequence`` describes either a protein or a DNA sequence.

    ``Sequence`` contains methods that act on single sequences, for example:

    *   tandem repeats


    Attributes:
        seq (str): The sequence.
        seq_standard_aa (str): The sequence with standard amino acids only
    """


    def __init__(self, seq, id = None):

        if not isinstance(seq, str):
                raise Exception('The seq value is not a String')
        self.seq = seq.upper()
        for i in self.seq:
            if i not in config['lAll_amino_acid']:
                raise Exception("{} is not in config['lAll_amino_acid']: {}".format(i, config['lAll_amino_acid']))

        self.seq_standard_aa = repeat.standardize(self.seq)

        if id:
            self.id = id

    def create(file, format):

        """ Create sequence(s) from file.

        Create sequence(s) from file.

        Args:
            file (str): Path to input file
            format (str):  Either "fasta" or "pickle"

        .. todo:: Write checks for ``format`` and ``file``.
        """

        if format == 'fasta':
            lSeq = sequence_io.read_fasta(file)
            return [Sequence(iSeq, iID) for iSeq, iID in lSeq]
        if format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        else:
            raise Exception("Output format {} is not implemented for sequence.write()".format(format))

    def write(self, file, format):

        """ Write sequence to file.

        Write sequence to file using one of two formats.

        Args:
            file (str): Path to output file
            format (str):  Either "fasta" or "pickle"

        .. todo:: Write checks for ``format`` and ``file``.

        """

        if format == 'fasta':
            sequence_io.write(self.seq, file)
        elif format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        else:
            raise Exception("Output format {} is not implemented for sequence.write()".format(format))

    def detect(self, lHMM = None, denovo = None, **kwargs):

        """ Detects tandem repeats on ``self.seq`` from 2 possible sources.

        A list of ``Repeat`` instances is created for tandem repeat detections on the
        sequence from two possible sources:

        * Sequence profile hidden Markov models ``HMM``
        * de novo detection algorithms.


        Args:
            hmm (HMM): A list of ``HMM`` instances.
            denovo (bool): boolean
            *args: Parameters fed to denovo TR prediction and/or Repeat instantiation.
                E.g. ``calc_score = True``

        Returns:
            A list of ``Repeat`` instances.

        .. todo:: Return ``repeat_list`` instance instead of a list of repeats

        """

        if lHMM:
            if not isinstance(lHMM, list):
                raise Exception('The lHMM value is not a list.')
            for iHMM in lHMM:
                if not isinstance(iHMM, hmm.HMM):
                    raise Exception('At least one list element in the lHMM value is '
                         'not a valid instance of the HMM class.')

            lRepeat = []
            for iHMM in lHMM:
                # Detect TRs on self.seq with hmm using the Viterbi algorithm.
                most_likely_path = iHMM.viterbi(self.seq)
                logging.debug(most_likely_path)
                unaligned_msa = hmm_viterbi.hmm_path_to_non_aligned_tandem_repeat_units(self.seq, most_likely_path, iHMM.lD)
                if len(unaligned_msa) > 1:
                    # Align the msa
                    aligned_msa = repeat_align.realign_repeat(unaligned_msa)
                    if len(aligned_msa) > 1:
                        # Create a Repeat() class with the new msa
                        if 'repeat' in kwargs:
                            lRepeat.append(repeat.Repeat(aligned_msa, **kwargs['repeat']))
                        else:
                            lRepeat.append(repeat.Repeat(aligned_msa))

            # Set begin coordinate for all repeats
            for iRepeat in lRepeat:
                self.repeat_in_sequence(iRepeat)

            return repeat_list.Repeat_list(lRepeat)

        elif lHMM == []:
            logging.debug("lHMM == []")
            return None

        elif denovo:
            if 'detection' in kwargs:
                lPredicted_repeat = repeat_detection_run.run_TRD([self], **kwargs['detection'])[0]
            else:
                lPredicted_repeat = repeat_detection_run.run_TRD([self])[0]

            lRepeat = []

            for jTRD,jlTR in lPredicted_repeat.items():
                for iTR in jlTR:
                    if 'repeat' in kwargs:
                        iTR = repeat.Repeat(iTR.msa, begin = iTR.begin, **kwargs['repeat'])
                    else:
                        iTR = repeat.Repeat(iTR.msa, begin = iTR.begin)


                    # Consider only tandem repeats that have a repeat unit predicted to be at least one character long.
                    if iTR.lD > 0:

                        # Save l, n, MSA, TRD, scores, sequence_type, position in sequence of given type
                        iTR.TRD = jTRD

                        # Sanity check repeat and set begin coordinate for all repeats
                        if not self.repeat_in_sequence(iTR):
                            logging.debug("The tandem repeat is not part of the sequence. Detector: {}".format(iTR.TRD))
                            continue

                        lRepeat.append(iTR)

            return repeat_list.Repeat_list(lRepeat)

        else:
            raise Exception("Either require denovo detection, or present an HMM")

    def set_repeat_list(self, repeat_list, tag):

        """ Add `repeat_list` as attribute to this `sequence` instance.

        Add `repeat_list` as attribute to this `sequence` instance. Access `repeat_list`
        as self.dRepeat_list[tag]

        Args:
            repeat_list (repeat_list): A repeat_list instance.
            tag (str): A identifier for the repeat_list
        """

        if not hasattr(self,"dRepeat_list"):
            self.dRepeat_list = {}

        self.dRepeat_list[tag] = repeat_list

    def annotate(self, data, tag):

        if not hasattr(self,"dAnnotations"):
            self.dAnnotations = {}

        self.dAnnotations[tag] = data


    def get_annotation(self, tag):

        if hasattr(self,"dAnnotations") and tag in self.dAnnotations:
            return self.dAnnotations[tag]
        else:
            return None

    def repeat_in_sequence(self, myRepeat):

        """ Sanity check whether the `repeat` is part of this `sequence`. In case,
        calculate the position of the `repeat` within the `sequence`.

        If yes: Return True, set repeat.begin to corrected value if necessary.
        If no: Return False.
        Perform sanity check on sequences where all amino acids are or are converted to
        standard amino acids.

        Args:
            sequence (sequence): A sequence instance.

        Returns:
            bool: True if repeat is part of sequence, else false

        .. todo:: Decide whether save_original_msa is needed here.
        """

        repeat_sequence = repeat.get_repeat_sequence(myRepeat.msa_standard_aa)
        starts = [m.start()+1 for m in re.finditer(repeat_sequence,self.seq_standard_aa)] # The first letter in the sequence is counted as 1 (not 0, as in python).

        if len(starts) != 0: # Is the tandem repeat predicted correctly?
            if not hasattr(myRepeat,"begin") or not myRepeat.begin in starts:
                myRepeat.begin = starts[0]
            myRepeat.save_original_msa(self.seq)
            return True
        else:
            return False
