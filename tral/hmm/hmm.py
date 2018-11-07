# (C) 2012-2015 Elke Schaper
# (C) 2014 Julia Pecerska

"""
    :synopsis: A cyclic hidden Markov model that describes sequence tandem repeats.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb.sib.ch>
"""

import datetime
import logging
import numpy as np
import os
import pickle
import shutil
import subprocess
import tempfile

from tral.hmm import hmm_io, hmm_viterbi
from tral.repeat.repeat import Repeat
from tral import configuration

LOG = logging.getLogger(__name__)

CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["hmm"]


################################### HMM class ############################
class HMM:

    """ Sequence profile hidden Markov models (HMMs) are used to model genomic
    sequence containing a tandem repeat.

    Note:
        The HMM implemented here is described in detail in

        Schaper, E., Gascuel, O. & Anisimova, M. Deep conservation of human
        protein tandem repeats within the eukaryotes.
        Molecular Biology and Evolution (2014).


        All notations used here are based on the notations introduced in

        Eddy, S. R. Profile hidden Markov models.
        Bioinformatics 14, 755–763 (1998).

    The HMM contains flanking states ("N", "C"), match states ("M1", "M2", ...)
    and insertion states ("I1", "I2", ..). As deletion states are non-emitting,
    they are merged with the other states.

    The HMM state "N" emits the sequence before the tandem repeat.
    The HMM state "C" emits the residual sequence after the tandem repeat.
    The HMM states "M1", "M2", ... can emit sequence within the tandem repeat.
    "M1" describes the first amino acid/nucleotide in the repeat unit,
    "M2" the second, and so on]. The HMM states "I1", "I2", can emit sequence
    representing insertions within the tandem repeat.

    The transition probabilities between states are stored in ``p_t``.
    The emission probabilities of every state are stored in ``p_e``.
    The probabilities for any state to be the null state are stored in ``p_0``.

    The structure and probabilities defining standard sequence profile HMMs can
    be

    * taken from HMM databases such as PFAM,
    * calculated from sequence alignments (here: an alignment of TR units) using for example the HMMER suite,
    * defined in house.

    All three ways are implemented here.

    The HMM is used to detect the maximum likelihood occurrence of a TR modeled
    by the HMM on a sequence. The Viterbi algorithm is implemented to do this
    detection.

    Attributes:
        id (str): The id describes the origin of the HMM.
                    E.g., PFAMID "PF00069".
        states (list of str): The names of all states that the HMM contains.
                              E.g. ["N", "M1", "I1", "M2", "I2", "C"]
        insertion_states (list of str): The names of all insertion states that
                                        the HMM contains. E.g. ["I1", "I2"]
        match_states (list of str): The names of all match states that the HMM
                                    contains. E.g. ["M1", "M2"]
        l_effective (int): The number of states in the HMM.
                           E.g. l_effective = 4
        alphabet (list of str): The letters that the HMM states can emit.
                                E.g. all amino acids in the LG matrix:
                                ["A", "C", ...]
        p_t (dict of dict of float): A dictionary
            {state0: {state1: transition probability, ...}, ...}
            between all `states` as log10.
        p_0 (dict of dict of float): A dictionary of the null probabilities to
                                    be in any of the `states`:
                                    {`states[0]`: 0.1, `states[1]`: 0.2, ...}
        p_e (dict of dict of float): A dictionary of emission probabilities for
                every state and every letter in `alphabet`:
                {`states[0]`: {`alphabet[0]`: emission probability, ...}, ...}
        hmmer (dict of dict of float): A dictionary of all information
                                       extracted from a Hmmer model.
    """

    def __str__(self):
        """
            Create string for HMM instance.
        """
        try:
            tmp = [self.p_e[i] for i in self.match_states]
            tmp = "".join([max(i.keys(), key=(lambda key: i[key]))
                           for i in tmp])
            hmm = "cpHMM ID: {}\ncpHMM length: {}\n".format(
                self.id,
                self.l_effective)
            hmm += "Most likely motif: {}".format(tmp)
        except (AttributeError, TypeError):
            hmm = "<HMM instance>"
            LOG.warning("Could not create string of HMM instance.")

        return hmm

    # Use functions defined in other methods as class functions.
    # Static function
    @staticmethod
    def read(hmm_filename, id=None):
        """Read HMM file in HMMER3 format.

        HMMER3 file format is described in detail in
        ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf,
        section 8.

        The general file format is::

            HMMER3/b [3.0b2 | June 2009]
            NAME fn3
            ACC PF00041.12
            DESC Fibronectin type III domain
            LENG 86
            ALPH amino
            (...)

            HMM          A       C       D       E       F       G       H       I    (...)    Y
                        m->m    m->i    m->d    i->m    i->i    d->m    d->d
              COMPO   2.70271 4.89246 3.02314 2.64362 3.59817 2.82566 3.74147 3.08574 (...) 3.22607
                      2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                      0.00338 6.08833 6.81068 0.61958 0.77255 0.00000       *
              1       3.16986 5.21447 4.52134 3.29953 4.34285 4.18764 4.30886 3.35801 (...) 3.93889 1 - -
                      2.68629 4.42236 2.77530 2.73088 3.46365 2.40512 3.72505 3.29365 (...) 3.61514
                      0.09796 2.38361 6.81068 0.10064 2.34607 0.48576 0.95510
              2       2.70230 5.97353 2.24744 2.62947 5.31433 2.60356 4.43584 4.79731 (...) 4.25623 3 - -
                      2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                      0.00338 6.08833 6.81068 0.61958 0.77255 0.48576 0.95510
              (...)
              86      3.03720 5.94099 3.75455 2.96917 5.26587 2.91682 3.66571 4.11840 (...) 4.99111 121 - E
                      2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                      0.00227 6.08723       * 0.61958 0.77255 0.00000       *
            //

        .. important:: The probabilities in this format are stored as negative
            natural log probabilities, e.g. -ln(0.25) = 1.38629. The special case
            of 0 probability is stored as ``*`` which in fact means -∞ (minus
            infinity).

        Args:
            hmm_filename (str): Path to the file with model data in the HMMER3
                file format.
            id (str, optional): The identifier for the model to be returned.
                E.g. for Pfam the ``id`` can look like this: 'PF00560'.
                If defined, the function returns only the HMM with an identifier
                that matches the provided id. If a model in the file does not have
                an identifier to check against, it is skipped.

        Returns:
            dict: A dictionary of parameters required to initialise the HMM
            including the id, alphabet, emission probabilities and transition
            probabilities.

            Output format::

                {
                    'id': 'PF08261.7',
                    'letters': ['A', 'C', 'D', ..., 'Y'],
                    'COMPO':
                        {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                          'emissions': [2.28205, 5.14899, ..., 1.92022],
                          'transitions': [0.01467, 4.62483, ..., -inf]},
                    '1':
                        {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                         'emissions': [1.00089, 4.54999, ..., 5.23581],
                         'transitions': [0.01467, 4.62483, ..., 0.95510]},
                    ...
                    '8':
                        {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                         'emissions': [4.12723, 5.39816, ..., 4.58094],
                         'transitions': [0.00990, 4.62006, ..., -inf]}
                }

        See Also:
            hmm_io.read()

        """
        return hmm_io.read(hmm_filename, id)

    # Function that take the HMM instance as argument:
    def viterbi(self, *args):
        return hmm_viterbi.viterbi(self, *args)

    def __init__(self, hmmer_probabilities, sequence_type="AA"):
        """ HMM class __init_ module.

        __init__ takes HMM parameters (including the alphabet, emission
        probabilities and transition probabilities) as input, and assigns
        them to class attributes.

        Args:
            hmmer_probabilities (dict): A dictionary with HMM parameters.
        """
        self.hmmer = hmmer_probabilities
        self.id = hmmer_probabilities['id']
        self.alphabet = self.hmmer['letters']
        self.sequence_type = sequence_type

        self.l_effective = max([int(key) for key
                               in hmmer_probabilities.keys() if key.isdigit()])

        # Initialise all HMM states to default value (e.g. transition to or
        # from terminal state all have the same cost 0).
        self.initialise_HMM_structure(self.l_effective)

        d_translate_states = {
            i: i[
                1:] for i in self.match_states +
            self.insertion_states}
        for i in self.terminal_states:
            d_translate_states[i] = "COMPO"

        # Store null model emission probabilities for ["N","C"]
        # from hmmer_probabilities["COMPO"]
        # Store match_state and insertion_state emission probabilities for
        # ["M0",...] from hmmer_probabilities["0"],...
        for i_state in self.match_states:
            i_emission_probabilities = self.hmmer[
                d_translate_states[i_state]]['emissions'][:len(self.alphabet)]
            self.set_emission_probability_hmmer3(
                i_state,
                i_emission_probabilities)
        for i_state in (self.insertion_states + self.terminal_states):
            i_emission_probabilities = \
                self.hmmer[d_translate_states[i_state]]['insertion_emissions'][:len(self.alphabet)]
            self.set_emission_probability_hmmer3(
                i_state,
                i_emission_probabilities)

        # Store transition probabilities
        # (m->m    m->i    m->d    i->m    i->i    d->m    d->d):
        # First: All state defined in HMMer model
        for i, i_state in enumerate(self.match_states):
            i_transition_probabilities = self.hmmer[
                d_translate_states[i_state]]['transition']
            # Completing the circle of the HMM.
            if i == self.l_effective - 1:
                i_transition_probabilities = self.set_circle_transition_probability_hmmer3(
                    i_transition_probabilities,
                    d_translate_states[i_state])
            LOG.debug("i_state: %s", i_state)
            LOG.debug(
                "set_circle_transition_probability_hmmer3: %s",
                i_transition_probabilities)
            self.set_transition_probability_hmmer3(
                i,
                i_transition_probabilities)

        LOG.debug("p_t before deletion: %s", self.p_t)

        # Translate deletion states into direct transitions from match to match
        # states.
        if self.l_effective > 1:
            p_t_d = {}
            for i, i_match_state in enumerate(self.match_states):
                p_t_d[i_match_state] = self.get_direct_transition_probabilities_for_deletions(
                    i,
                    i_match_state)
            for i_match_state in self.match_states:
                for i_goal_state, iP in p_t_d[i_match_state].items():
                    self.p_t[i_match_state][i_goal_state] = iP

        LOG.debug("p_t after deletion: %s", self.p_t)

        # As it is never likely to got from a terminal state to a deletion
        #  tate or vice versa, this transition probabilities can be ignored.
        # Thus, you may advance and:
        # Delete everything you ever knew about deletion states.

        self.states = [self.terminal_states[
            0]] + self.match_states + self.insertion_states + [self.terminal_states[1]]
        self.p_t = {i_state: {i: j for i, j in self.p_t[
            i_state].items() if i not in self.deletion_states} for i_state in self.states}
        self.p_0 = {
            i: j for i,
            j in self.p_0.items() if i not in self.deletion_states}

        del self.deletion_states

    @staticmethod
    def create(input_format, file=None, repeat=None, **kwargs):
        """ Creates a HMM instance from 2 possible input formats.

        A `HMM` instance is created from one of the two possible inputs:

        * HMMER3 model
        * ``Repeat`` instance

        Args:
            input_format (str): The file format of the input file
            file (str): Path to the file containing the HMM parameters
                encoded in the HMMER3 format.
            repeat (Repeat): A Repeat object with an MSA that can be
                transformed into an HMM.

        Returns:
            `HMM`: An initialized instance of the `HMM` class.

        Raises:
            Exception: if the HMMER3 file provided does not exist.
            Exception: if the repeat value provided is not an instance of the
                Repeat class.
            Exception: if no parameters are provided to create the HMM from.

        .. warning:: Will initialize all HMM's in a file but will only return
            the first one.
        .. todo:: Fix the previous warning (agree on the way how it should be
            done first)

        """

        if input_format == 'hmmer':
            # Read only first element from HMMER3 file.
            if not os.path.exists(file):
                raise Exception('HMMER3 file does not exist.')
            hmmer_probabilities = next(HMM.read(file))
            sequence_type = "AA"
        elif input_format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        elif input_format == 'repeat':
            if not isinstance(repeat, Repeat):
                raise Exception('The repeat value is not a valid instance of '
                                'the Repeat class.')
            if 'hmmbuild' in kwargs:
                hmmer_probabilities = HMM.create_from_repeat(repeat, **kwargs['hmmbuild'])
            else:
                hmmer_probabilities = HMM.create_from_repeat(repeat)
            sequence_type = repeat.sequence_type

        else:
            raise Exception("Unknown input format: {}.".format(input_format))

        LOG.debug(hmmer_probabilities)
        return HMM(hmmer_probabilities, sequence_type)

    @staticmethod
    def create_from_repeat(tandem_repeat, hmm_copy_path=None,
                           hmm_copy_id=None):
        """ Get HMM parameters (including the alphabet, emission probabilities
        and transition probabilities) from a tandem repeat.

        An HMM is created from ``tandem_repeat`` using :command:`hmmbuild`
        and is saven in the HMMER3 file format.
        Next, HMM parameters are retrieved from the HMMER3 file and returned.

        Args:
            tandem_repeat (TR):  A Repeat class instance of the TR to be
                transformed into an HMM.
            hmm_copy_path (str): Path to where a copy of the created HMM
                will be stored.
                If None, no copies are saved.
            hmm_copy_id (str): HMM id which will serve as the filename in case
                a copy should be saved.
                If None and ``hmm_copy_path`` is specified the name is
                generated from current date and time combination.

        Returns:
            HMM parameters (dict): For a sample structure please refer to
            :py:hmm_io.read:.

        .. todo:: Make sure the id is correctly read in.

        """
        timestamp = datetime.datetime.now()
        tmp_id = timestamp.strftime('%Y:%m:%d:%H:%M:%S:%f')

        # Create a temporary directory
        tmp_dir = tempfile.mkdtemp()

        # Save TR as Stockholm file
        stockholm_file = os.path.join(tmp_dir, tmp_id + ".sto")

        # O repeat_io.save_repeat_stockholm(tandem_repeat.msaD, stockholm_file)
        tandem_repeat.write(file=stockholm_file, file_format="stockholm")

        # Define hmmbuild sequence type flag.
        if tandem_repeat.sequence_type == "AA":
            sequence_type_flag = "--amino"
        elif tandem_repeat.sequence_type == "RNA":
            sequence_type_flag = "--rna"
        else:
            sequence_type_flag = "--dna"

        # Run HMMbuild to build a HMM model, and read model
        p = subprocess.Popen([CONFIG["hmmbuild"], sequence_type_flag, tmp_id + ".hmm",
                              tmp_id + ".sto"],
                             stdout=subprocess.PIPE, stderr=None, cwd=tmp_dir)
        p.wait()

        hmm_file = os.path.join(tmp_dir, tmp_id + ".hmm")
        if hmm_copy_path:
            if not os.path.exists(hmm_copy_path):
                LOG.critical('Specified path for the file copy does not '
                             'exist, not saving copy: %s.', hmm_copy_path)
            else:
                if hmm_copy_id:
                    shutil.copy(hmm_file, os.path.join(hmm_copy_path,
                                                       hmm_copy_id + ".hmm"))
                else:
                    shutil.copy(hmm_file, hmm_copy_path)

        hmmer_probabilities = next(
            HMM.read(
                hmm_filename=hmm_file,
                id=hmm_copy_id))

        shutil.rmtree(tmp_dir)

        return hmmer_probabilities

    def write(self, file, output_format, *args):
        """ Write ``HMM`` to file.

        Write ``HMM`` to file. Currently, only pickle is implemented.

        Args:
            file (str): Path to input file
            output_format (str):  Either "fasta", "pickle" or "stockholm"

        .. todo:: Write checks for ``format`` and ``file``.

        """

        if output_format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        else:
            raise Exception('output_format is unknown: {}'.format(output_format))

    def hmm_example(self):
        """ .. todo:: Is this method needed?
        """
        # states = ["N", "B", "M1", "M2", "M3", "E", "C"]
        self.states = ["N", "M1", "M2", "C"]

        # Initialisation
        self.p_t = {iS: {iS2: 0 for iS2 in self.states}
                    for iS in self.states[:-1]}
        self.p_0 = {iS: 1 / 3 for iS in self.states}

        # Transition probabilities
        # Feed Values to p_t
        self.p_t["N"] = {"N": 0.5, "M1": 0.5}
        self.p_t["M1"] = {"M2": 0.5, "C": 0.5}
        self.p_t["M2"] = {"M1": 0.5, "C": 0.5}
        self.p_t["C"] = {"C": 1}

        # emissions
        self.emissions = ["A", "C", "G", "T"]

        # emission probabilities
        self.p_e = {iS: {iE: 0.25 for iE in self.emissions}
                    for iS in self.states}
        self.p_e['M1'] = {"A": 0.9, "C": 0.025, "G": 0.025, "T": 0.025}
        self.p_e['M2'] = {"A": 0.025, "C": 0.9, "G": 0.025, "T": 0.025}

    def initialise_HMM_structure(self, l_effective):
        """ Initialise the HMM with all states and standard values for
            transition and emission probabilities.

        Initialise all states
        Set transition probabilities to None/0
        Set initial probability for all states equal.

        Args:
          l_effective (int): The number of match states in the HMM.

        """

        # Build a HMM from these likelihoods.
        self.match_states = ["M{0}".format(str(i + 1)) for i in range(l_effective)]
        self.insertion_states = ["I{0}".format(str(i + 1)) for i in range(l_effective)]
        self.deletion_states = ["D{0}".format(str(i + 1)) for i in range(l_effective)]
        self.terminal_states = ["N", "C"]
        self.states = [self.terminal_states[0]] + self.match_states + \
            self.insertion_states + self.deletion_states + [self.terminal_states[1]]
        LOG.debug("HMM states: %s", self.states)

        # Initialisation
        # The transition probability is initially set to None for all states.
        self.p_t = {iS: {iS2: None for iS2 in self.states}
                    for iS in self.states}
        # The initial probability is equal in all states. (np.log10(1)=0)
        self.p_0 = {iS: 0 for iS in self.states}

        # Transition probabilities

        # First, mark which transitions are not penalised,
        # that is bear a 'probability of np.log10(1)=0':
        # - Transitions from "N" to any "M_k"
        # - Transitions from any "M_k" to "C"
        # - Transitions from "M_k" to "M_k+1" (for k`n) and "M_k" to "M_0" (for k==n)
        # - Transitions from "N" to "N"
        # - Transitions from "C" to "C"
        for i in range(l_effective):
            self.p_t["N"]["M{0}".format(i % l_effective + 1)] = 0
            self.p_t["M{0}".format(i % l_effective + 1)]["C"] = 0
            self.p_t["M{0}".format(i % l_effective + 1)]["M{0}".format((i + 1) % l_effective + 1)] = 0
        self.p_t["N"]["N"] = 0
        self.p_t["C"]["C"] = 0

        self.p_e = {}

    def set_emission_probability_hmmer3(self, state, l_emission_probabilities):
        """ Convert and set emission probabilities from HMMER3 file to class
            attribute.

        Set p_e of states to `l_emission_probabilities` for `state` given
        `self.alphabet`.
        In HMMER3 data, emission probabilities are -ln(p).

        It is e.g.::

            self.alphabet =
            ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
             'R', 'S', 'T', 'V', 'W', 'Y']

        Return log10(p), that is convert between the two. Conversion:
        p_Local = - p_HMM * log10(e)

        Args:
            state (str): The state of the HMM. E.g. the fourth match state 'M4'.
            l_emission_probabilities (list of float): E.g.
                ['3.27687', '2.31397', '3.32252', '3.12746', '2.89175',
                '3.34719', '2.28730', '3.54139', '2.53154', '2.64774',
                '3.75733', '3.23860', '3.57894', '3.26290', '2.66343',
                '2.61544', '2.91770', '3.26739', '5.09378', '3.16816']

        """

        l_emission_probabilities = \
            [-(i) * np.log10(np.exp(1)) for i in l_emission_probabilities]
        self.p_e[state] = \
            {iL: iEP for iL, iEP in
             zip(self.alphabet, l_emission_probabilities)}

    def set_circle_transition_probability_hmmer3(
            self,
            l_transition_probabilities,
            final_state):
        """ Calculate transition probabilites between states that exist due to
            the circularity of the tandem repeat HMM.

        In HMMER3 models, no transition in circle from the last match state to
        the first match state are included.

        Values need to be found for the following transitions:
        m->m    m->i    m->d    i->m    i->i    d->m    d->d

        Compare ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf
        pp.84-85.

        Input: -ln(p) Output: Also -ln(p')

        When the HMM only has a single match state, i.e. l_effective == 1, the
        final match state equals the only match state, and no transition
        probabilities can be calculated via averaging over all other states.
        Therefore, the its transition probabilities remain unchanged.

        Args:
          l_transition_probabilities (list of float): The first parameter.
          final_state (str): The second parameter.

        Returns:
          list of float: transition probabilities

        """

        if self.l_effective == 1:
            return l_transition_probabilities

        # In self.hmmer, besides the probabilities for states 1, 2, 3, ... data
        # for "COMPO", 'letters', and 'id' are saved and should be ignored in
        # the following.
        # Also, we ignore the transition probabilities for final_state, as our
        # goal here is to recalculate these based on the averages of the
        # transition probabilities from all other states.
        skip_states = ['COMPO', 'letters', final_state, 'id']

        # M->M, M->I, M->D: Use an average from all other transitions from
        # match states.
        LOG.debug("self.hmmer.keys(): %s", self.hmmer.items())
        LOG.debug("self.hmmer.items(): %s", self.hmmer.items())
        mean_M_M = np.mean([np.exp(-iP["transition"][0])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        mean_M_I = np.mean([np.exp(-iP["transition"][1])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        mean_M_D = np.mean([np.exp(-iP["transition"][2])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])

        sum_p = np.sum([mean_M_M, mean_M_I, mean_M_D])
        l_transition_probabilities[
            :3] = [-np.log(i / sum_p) for i in [mean_M_M, mean_M_I, mean_M_D]]

        # I->M and I->I: Use an average from all other transitions from
        # insertion states.
        mean_I_M = np.mean([np.exp(-iP["transition"][3])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        mean_I_I = np.mean([np.exp(-iP["transition"][4])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        sum_p = np.sum([mean_I_M, mean_I_I])
        l_transition_probabilities[
            3:5] = [-np.log(i / sum_p) for i in [mean_I_M, mean_I_I]]

        # D->M and D->D: Use an average from all other transitions from
        # deletion states.
        mean_D_M = np.mean([np.exp(-iP["transition"][5])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        mean_D_D = np.mean([np.exp(-iP["transition"][6])
                            for i_state, iP in self.hmmer.items() if i_state not in skip_states])
        sum_p = np.sum([mean_D_M, mean_D_D])
        l_transition_probabilities[
            5:] = [-np.log(i / sum_p) for i in [mean_D_M, mean_D_D]]

        return l_transition_probabilities

    def set_transition_probability_hmmer3(
            self,
            state_index,
            l_transition_probabilities):
        """ Convert and assign transition probabilities from a HMMER3 file to
            HMM attribute `p_t`.


        Convert and assign transition probabilities from a HMMER3 file to
        HMM attribute `p_t`.
        Set ``p_t`` of states to ``l_transition_probabilities`` for the
        ``state_index`` th state.

        In HMMER3 data, transition probabilities are -ln(p).
        Return log10(p), i.d. convert between the two.
        Conversion: p_Local = - p_HMM * log10(e)

        Args:
            state_index (int): E.g. 4
            l_transition_probabilities (list of float): A list of transition
                probabilities as found in a single row for a certain state in
                HMMER3 files. E.g. ['0.00021', '8.85487', '9.57722', '0.61958',
                '0.77255', '0.48576', '0.95510']
        """

        l_transition_probabilities = [-(i) * np.log10(np.exp(1)) for i
                                      in l_transition_probabilities]

        # Match state state_index -> Match state state_index + 1
        iM = self.match_states[state_index]
        iM_1 = self.match_states[(state_index + 1) % self.l_effective]
        self.p_t[iM][iM_1] = l_transition_probabilities[0]

        # Match state state_index -> Insertion state state_index
        self.p_t[iM][
            self.insertion_states[state_index]] = l_transition_probabilities[1]

        # Match state state_index -> Deletion state state_index + 1
        iD_1 = self.deletion_states[(state_index + 1) % self.l_effective]
        self.p_t[iM][iD_1] = l_transition_probabilities[2]

        # Insertion state state_index ->  Match state state_index + 1
        iI = self.insertion_states[state_index]
        self.p_t[iI][iM_1] = l_transition_probabilities[3]

        # Insertion state state_index ->  Insertion state state_index
        iI = self.insertion_states[state_index]
        self.p_t[iI][iI] = l_transition_probabilities[4]

        # Deletion state state_index -> Match state state_index + 1
        iD = self.deletion_states[state_index]
        self.p_t[iD][iM_1] = l_transition_probabilities[5]

        # Deletion state state_index -> Deletion state state_index + 1
        self.p_t[iD][iD_1] = l_transition_probabilities[6]

    def get_direct_transition_probabilities_for_deletions(
            self,
            state_index,
            state):
        """ Calculate updated transition probabilities for HMMs without
            deletion states.

        Calculate updated transition probabilities for HMMs without deletion
        states. Deletion states are non-emitting states. Therefore, they can be
        merged into the remaining emitting states. E.g. M -> D -> I would be
        equivalent to M -> I. The total transition probability from M to I
        would then be the product of the transition probabilities from M to D
        and from D to I (equivalent to the sum of the log10 probabilities.).

        This method translates all deletions following `state` with
        `state_index` into direct match state to match state transitions, and
        returns the new transition probabilities.

        Args:
            state_index (int): The first parameter.
            state (str): The second parameter.

        Returns:
            {state1: float}: The transition probability from `state` with
            `state_index` to `state1` is calculated for all (?) possible
            `state1`.
        """

        transition = {}
        p_deletion_state = self.p_t[state][
            self.deletion_states[
                (state_index + 1) %
                self.l_effective]]
        for i_state_index in range(1, self.l_effective):
            i_deleted_state_index = (state_index + i_state_index) % self.l_effective
            i_deleted_state = self.deletion_states[i_deleted_state_index]
            i_goal_state_index = (state_index + i_state_index + 1) % self.l_effective
            i_goal_state = self.match_states[i_goal_state_index]
            transition[i_goal_state] = p_deletion_state + \
                self.p_t[i_deleted_state][i_goal_state]
            p_deletion_state += self.p_t[i_deleted_state][
                self.deletion_states[i_goal_state_index]]

        return transition


############################### External parameters ######################

def hmmer3_emission_probabilities(hmmer_probabilities, letters, l_match):
    '''
    Get emission probabilities from hmmer3 hmm file.
    In hmm file, emission probabilities are -ln(p).
    Return log10(p), i.d. convert between the two.
    Conversion: p_Local = - p_HMM * log10(e)

    Parameters (e.g.)::

        letters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N',
                   'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
        l_match = ['M'+str(i) for i in range(24)]

    Return format (pseudo code):
    [{iA: np.log10(p(iA,iM)) for iA in alphabet.keys()} for iM in l_match]

    '''

    # Test: Is the number of match states in both models equal?
    if not len(hmmer_probabilities.keys()) - 2 == len(l_match):
        print('Match states HMMER: {0} Match states local: {1}'.format(
            len(hmmer_probabilities.keys()) - 2, len(l_match)))
        raise ValueError('The number of match states in HMMer model and local'
                         'model does not match')

    # Test: Are all `letters` represented in the HMMER HMM?
    if any(iL not in hmmer_probabilities['letters'] for iL in letters):
        missing = [
            iL for iL in letters if iL not in hmmer_probabilities['letters']]
        print('Missing representation in Hmmer File: {0}'.format(missing))
        raise ValueError('Some letters in the local HMM are not represented',
                         'in the HMMER HMM.')

    return [{iL: -(iP) * np.log10(np.exp(1)) for iL,
             iP in zip(hmmer_probabilities['letters'],
                       data['emissions']) if iL in letters} for key,
            data in hmmer_probabilities.items() if key not in ['letters',
                                                               'COMPO']]
