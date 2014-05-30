
from tandemrepeats.repeat import repeat, repeat_realign
from tandemrepeats.hmm import hmm, hmm_viterbi



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

    def detect(self, lHMM = None, denovo = None, *args):

        """ Detects tandem repeats on ``self.seq`` from 2 possible sources.

        A list of ``Repeat`` instances is created for tandem repeat detections on the
        sequence from two possible sources:

        * Sequence profile hidden Markov models ``HMM``
        * de novo detection algorithms.


        Args:
            hmm (HMM): A list of ``HMM`` instances.
            denovo (list of str): A list of tandem repeat detection algorithm names,
                e.g. ``["t-reks", "xstream"]``
            *args: Parameters fed to the Repeat instantiation. E.g. ``calc_score = True``

        Returns:
            A list of ``Repeat`` instances.

        .. todo:: Return ``repeat_list`` instance instead of a list of repeats

        """

        if lHMM:
            if not isinstance(lHMM, list):
                raise Exception('The lHMM value is not a list.')
            for iHMM in lHMM:
                if not isinstance(iHMM, HMM):
                    raise Exception('At least one list element in the lHMM value is '
                         'not a valid instance of the HMM class.')

            lRepeat = []
            for i in iHMM:
                # Detect TRs on self.seq with hmm using the Viterbi algorithm.
                most_likely_path = hmm.viterbi(self.seq)
                unaligned_msa = hmm_viterbi.hmm_path_to_non_aligned_tandem_repeat_units(self.seq, most_likely_path, hmm.lD)
                if len(unaligned_msa) > 1:
                    # Align the msa
                    aligned_msa = repeat_realign.realign_repeat(unaligned_msa)
                    if len(aligned_msa) > 1:
                        # Create a Repeat() class with the new msa
                        lRepeat.append(repeat.Repeat(aligned_msa, *args))
            return lRepeat

        elif denovo:
            raise Exception('Not implemented yet.')
        else:
            raise Exception("Neither of the required inputs provided!")
