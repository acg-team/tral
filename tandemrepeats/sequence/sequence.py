import logging

from tandemrepeats.repeat import repeat, repeat_align
from tandemrepeats.repeat_list import repeat_list
from tandemrepeats.hmm import hmm, hmm_viterbi
from tandemrepeats.sequence import repeat_detection_run


logging.basicConfig()
logger = logging.getLogger(__name__)

class Sequence:

    """ A ``Sequence`` describes either a protein or a DNA sequence.

    ``Sequence`` contains methods that act on single sequences, for example:

    *   tandem repeats


    Attributes:
        seq (str): The sequence.
    """


    def __init__(self, seq):

        if not isinstance(seq, str):
                raise Exception('The seq value is not a String')
        self.seq = seq

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
                        lRepeat.append(repeat.Repeat(aligned_msa, *args))
            return repeat_list.Repeat_list(lRepeat)

        elif denovo:

            if 'detection' in kwargs:
                lPredicted_repeat = repeat_detection_run.run_TRD([self.seq], **kwargs['detection'])[0]
            else:
                lPredicted_repeat = repeat_detection_run.run_TRD([self.seq])[0]

            lRepeat = []

            for jTRD,jlTR in lPredicted_repeat.items():
                for iTR in jlTR:
                    if 'repeat' in kwargs:
                        iTR = repeat_info.Repeat(iTR.msa, begin = iTR.begin, **kwargs['repeat'])
                    else:
                        iTR = repeat_info.Repeat(iTR.msa, begin = iTR.begin)


                    # Consider only tandem repeats that have a repeat unit predicted to be at least one character long.
                    if iTR.lD > 0:

                        # Save l, n, MSA, TRD, scores, sequence_type, position in sequence of given type
                        iTR.TRD = jTRD

                        repeat_in_sequence = iTR.repeat_in_sequence(self.seq, save_original_msa = True)
                        if not repeat_in_sequence:
                            logging.debug("The tandem repeat is not part of the sequence. Detector: {}".format(iTR.TRD))
                            continue

                        lRepeat.append(iTR)

            return repeat_list.Repeat_list(lRepeat)

        else:
            raise Exception("Either require denovo detection, or present an HMM")
