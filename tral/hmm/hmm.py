# (C) 2012-2015 Elke Schaper
# (C) 2014 Julia Pecerska

"""

    :synopsis: A cyclic hidden Markov model that describes sequence tandem repeats.

    .. moduleauthor:: Elke Schaper <elke@inf.ethz.ch>

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
from tral.repeat import repeat_io
from tral.repeat.repeat import Repeat
from tral import configuration

log = logging.getLogger(__name__)

c = configuration.Configuration.Instance()
config_general = c.config
config = config_general["hmm"]

################################### HMM class ############################


class HMM:

    """ Sequence profile hidden Markov models (HMMs) are used to model genomic
    sequence containing a tandem repeat.

    Note:
        The HMM implemented here is described in detail in

        Schaper, E., Gascuel, O. & Anisimova, M. Deep conservation of human protein tandem
        repeats within the eukaryotes. Molecular Biology and Evolution (2014).


        All notations used here are based on the notations introduced in

        Eddy, S. R. Profile hidden Markov models. Bioinformatics 14, 755–763 (1998).

    The HMM contains flanking states ("N",  "C"), match states ("M1", "M2", ...) and
    insertion states ("I1", "I2", ..). As deletion states are non-emitting, they are
    merged with the other states.

    The HMM state "N" emits the sequence before the tandem repeat.
    The HMM state "C" emits the residual sequence after the tandem repeat.
    The HMM states "M1", "M2", ... can emit sequence within the tandem repeat. "M1"
    describes the first amino acid/nucleotide in the repeat unit, "M2" the second, and so on.
    The HMM states "I1", "I2", can emit sequence representing insertions within the tandem repeat.

    The transition probabilities between states are stored in ``p_t``.
    The emission probabilities of every state are stored in ``p_e``.
    The probabilities for any state to be the null state are stored in ``p_0``.

    The structure and probabilities defining standard sequence profile HMMs can be

    * taken from HMM databases such as PFAM,
    * calculated from sequence alignments (here: an alignment of TR units) using for example the HMMER suite,
    * defined in house.

    All three ways are implemented here.

    The HMM is used to detect the maximum likelihood occurrence of a TR modeled by the HMM on a sequence.
    The Viterbi algorithm is implemented to do this detection.

    Attributes:
        id (str): The id describes the origin of the HMM. For example, if may be the PFAMID "PF00069".
        states (list of str): The names of all states that the HMM contains. E.g. ["N", "M1", "I1", "M2", "I2", "C"]
        insertion_states (list of str): The names of all insertion states that the HMM contains. E.g. ["I1", "I2"]
        match_states (list of str): The names of all match states that the HMM contains. E.g. ["M1", "M2"]
        lD (int): The number of states in the HMM. E.g. lD = 4
        alphabet (list of str): The letters that the HMM states can emit. E.g. all amino acids in the LG matrix. ["A", "C", ...]
        p_t (dict of dict of float): A dictionary {state0: {state1: transition probability, ...}, ...} between all `states` as log10.
        p_0 (dict of dict of float): A dictionary of the null probabilities to be in any of the `states`: {`states[0]`: 0.1, `states[1]`: 0.2, ...}
        p_e (dict of dict of float): A dictionary of emission probabilities for every state and every letter in `alphabet`: {`states[0]`: {`alphabet[0]`: emission probability, ...}, ...}
        hmmer (dict of dict of float): A dictionary of all information extracted from a Hmmer model.
    """

    def __str__(self):
        """
            Create string for HMM instance.
        """
        try:
            tmp = [self.p_e[i] for i in self.match_states]
            tmp = "".join([max(i.keys(), key=(lambda key: i[key]))
                           for i in tmp])
            hmm = "cpHMM ID:          {}\ncpHMM length: {}\n".format(
                self.id,
                self.lD)
            hmm += "Most likely motif:      {}".format(tmp)
        except:
            hmm = "<HMM instance>"
            log.warning("Could not create string of HMM instance.")

        return hmm

    # Use functions defined in other methods as class functions.
    # Static function (declare as @static ?):
    read = hmm_io.read

    # Function that take the HMM instance as argument:
    def viterbi(self, *args):
        return hmm_viterbi.viterbi(self, *args)

    def __init__(self, hmmer_probabilities):
        """ HMM class __init_ module.

        __init__ takes HMM parameters (including the alphabet, emission probabilities
        and transition probabilities) as input, and assigns them to class attributes.

        Args:
            hmmer_probabilities (dict): A dictionary with HMM parameters.
        """
        self.hmmer = hmmer_probabilities
        self.id = hmmer_probabilities['id']
        self.alphabet = self.hmmer['letters']

        self.lD = max([int(key) for key in hmmer_probabilities.keys()
                       if key.isdigit()])

        # Initialise all HMM states to default value (e.g. transition to or
        # from terminal state all have the same cost 0).
        self.initialise_HMM_structure(self.lD)

        dTranslate_States = {
            i: i[
                1:] for i in self.match_states +
            self.insertion_states}
        for i in self.terminal_states:
            dTranslate_States[i] = "COMPO"

        # Store null model emission probabilities for ["N","C"] from hmmer_probabilities["COMPO"]
        # Store match_state and insertion_state emission probabilities for
        # ["M0",...] from hmmer_probabilities["0"],...
        for iState in self.match_states:
            iEmission_Probabilities = self.hmmer[
                dTranslate_States[iState]]['emissions'][
                :len(
                    self.alphabet)]
            self.set_emission_probability_hmmer3(
                iState,
                iEmission_Probabilities)
        for iState in (self.insertion_states + self.terminal_states):
            iEmission_Probabilities = self.hmmer[
                dTranslate_States[iState]]['insertion_emissions'][
                :len(
                    self.alphabet)]
            self.set_emission_probability_hmmer3(
                iState,
                iEmission_Probabilities)

        # Store transition probabilities (m->m    m->i    m->d    i->m    i->i    d->m    d->d):
        # First: All state defined in HMMer model
        for i, iState in enumerate(self.match_states):
            iTransition_Probabilities = self.hmmer[
                dTranslate_States[iState]]['transition']
            # Completing the circle of the HMM.
            if i == self.lD - 1:
                iTransition_Probabilities = self.set_circle_transition_probability_hmmer3(
                    iTransition_Probabilities,
                    dTranslate_States[iState])
            log.debug("iState: %s", iState)
            log.debug(
                "set_circle_transition_probability_hmmer3: %s",
                iTransition_Probabilities)
            self.set_transition_probability_hmmer3(
                i,
                iTransition_Probabilities)

        log.debug("p_t before deletion: %s", self.p_t)

        # Translate deletion states into direct transitions from match to match
        # states.
        if self.lD > 1:
            p_t_d = {}
            for i, iMatch_state in enumerate(self.match_states):
                p_t_d[iMatch_state] = self.get_direct_transition_probabilities_for_deletions(
                    i,
                    iMatch_state)
            for iMatch_state in self.match_states:
                for iGoal_state, iP in p_t_d[iMatch_state].items():
                    self.p_t[iMatch_state][iGoal_state] = iP

        log.debug("p_t after deletion: %s", self.p_t)

        # As it is never likely to got from a terminal state to a deletion state or vice versa, this transition probabilities can be ignored.
        # Thus, you may advance and:
        # Delete everything you ever knew about deletion states.

        self.states = [self.terminal_states[
            0]] + self.match_states + self.insertion_states + [self.terminal_states[1]]
        self.p_t = {iState: {i: j for i, j in self.p_t[
            iState].items() if not i in self.deletion_states} for iState in self.states}
        self.p_0 = {
            i: j for i,
            j in self.p_0.items() if i not in self.deletion_states}

        del self.deletion_states

    def create(format, file=None, repeat=None):
        """ Creates a HMM instance from 2 possible input formats.

        A `HMM` instance is created from one of the two possible inputs:

        * HMMER3 model
        * ``Repeat`` instance

        Args:
            hmmer_file (str): Path to the file containing the HMM parameters
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

        if format == 'hmmer':
            # Read only first element from HMMER3 file.
            if not os.path.exists(file):
                raise Exception('HMMER3 file does not exist.')
            hmmer_probabilities = next(HMM.read(file))
        elif format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        elif format == 'repeat':
            if not isinstance(repeat, Repeat):
                raise Exception('The repeat value is not a valid instance of '
                                'the Repeat class.')
            hmmer_probabilities = HMM.create_from_repeat(repeat)
        else:
            raise Exception("Neither of the required inputs provided!")

        log.debug(hmmer_probabilities)
        return HMM(hmmer_probabilities)

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

        .. todo:: Julia: Has to be tested and fixed!
        .. todo:: Julia: Resolve the id question for read.

        """
        timestamp = datetime.datetime.now()
        tmp_id = timestamp.strftime('%Y:%m:%d:%H:%M:%S:%f')

        # Create a temporary directory
        tmp_dir = tempfile.mkdtemp()

        # Save TR as Stockholm file
        stockholm_file = os.path.join(tmp_dir, tmp_id + ".sto")

        # O repeat_io.save_repeat_stockholm(tandem_repeat.msaD, stockholm_file)
        tandem_repeat.write(file=stockholm_file, format="stockholm")

        # Run HMMbuild to build a HMM model, and read model
        p = subprocess.Popen([config["hmmbuild"], "--amino", tmp_id + ".hmm",
                              tmp_id + ".sto"],
                             stdout=subprocess.PIPE, stderr=None, cwd=tmp_dir)
        p.wait()

        hmm_file = os.path.join(tmp_dir, tmp_id + ".hmm")
        if hmm_copy_path:
            if not os.path.exists(hmm_copy_path):
                log.critical('Specified path for the file copy does not '
                             'exist, not saving copy: %s.', hmm_copy_path)
            else:
                if hmm_copy_id:
                    shutil.copy(hmm_file, os.path.join(hmm_file_copy,
                                                       hmm_copy_id + ".hmm"))
                else:
                    shutil.copy(hmm_file, hmm_copy_path)

        hmmer_probabilities = next(
            HMM.read(
                hmm_filename=hmm_file,
                id=hmm_copy_id))

        shutil.rmtree(tmp_dir)

        return hmmer_probabilities

    def write(self, file, format, *args):
        """ Write ``HMM`` to file.

        Write ``HMM`` to file. Currently, only pickle is implemented.

        Args:
            file (str): Path to input file
            format (str):  Either "fasta", "pickle" or "stockholm"

        .. todo:: Write checks for ``format`` and ``file``.

        """

        if format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        else:
            raise Exception('format is unknown.')

    def HMM_example(self):
        """ Hi Jewels! Do we want such a function?
        """
        #states = ["N", "B", "M1", "M2", "M3", "E", "C"]
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

    def initialise_HMM_structure(self, lD):
        """ Initialise the HMM with all states and standard values for transition and
            emission probabilities.

        Initialise all states
        Set transition probabilities to None/0
        Set initial probability for all states equal.

        Args:
          lD (int): The number of match states in the HMM.

        """

        # Build a HMM from these likelihoods.
        self.match_states = ["M{0}".format(str(i + 1)) for i in range(lD)]
        self.insertion_states = ["I{0}".format(str(i + 1)) for i in range(lD)]
        self.deletion_states = ["D{0}".format(str(i + 1)) for i in range(lD)]
        self.terminal_states = ["N", "C"]
        self.states = [self.terminal_states[0]] + self.match_states + \
            self.insertion_states + self.deletion_states + [self.terminal_states[1]]
        log.debug("HMM states: %s", self.states)

        # Initialisation
        # The transition probability is initially set to None for all states.
        self.p_t = {iS: {iS2: None for iS2 in self.states}
                    for iS in self.states}
        # The initial probability is equal in all states. (np.log10(1)=0)
        self.p_0 = {iS: 0 for iS in self.states}

        # Transition probabilities

        # First, mark which transitions are not penalised, e.g. bear a 'probability of np.log10(1)=0':
        # - Transitions from "N" to any "M_k"
        # - Transitions from any "M_k" to "C"
        # - Transitions from "M_k" to "M_k+1" (for k`n) and "M_k" to "M_0" (for k==n)
        # - Transitions from "N" to "N"
        # - Transitions from "C" to "C"
        for i in range(lD):
            self.p_t["N"]["M{0}".format(i % lD + 1)] = 0
            self.p_t["M{0}".format(i % lD + 1)]["C"] = 0
            self.p_t["M{0}".format(i %
                                   lD + 1)]["M{0}".format((i + 1) %
                                                          lD + 1)] = 0
        self.p_t["N"]["N"] = 0
        self.p_t["C"]["C"] = 0

        self.p_e = {}

    def set_emission_probability_hmmer3(self, state, lEmission_Probabilities):
        """ Convert and set emission probabilities from HMMER3 file to class attribute.

        Set p_e of states to `lEmission_Probabilities` for `state` given `self.alphabet`
        In HMMER3 data, emission probabilities are -ln(p).

        It is e.g.::
            self.alphabet =
            ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        Return log10(p), that is convert between the two. Conversion:
        p_Local = - p_HMM * log10(e)

        Args:
            state (str): The state of the HMM. E.g. the fourth match state 'M4'.
            lEmission_Probabilities (list of float): E.g.
                ['3.27687', '2.31397', '3.32252', '3.12746', '2.89175', '3.34719', '2.28730', '3.54139', '2.53154', '2.64774', '3.75733', '3.23860', '3.57894', '3.26290', '2.66343', '2.61544', '2.91770', '3.26739', '5.09378', '3.16816']

        """

        lEmission_Probabilities = [-(i) * np.log10(np.exp(1))
                                   for i in lEmission_Probabilities]
        self.p_e[state] = {
            iL: iEP for iL,
            iEP in zip(
                self.alphabet,
                lEmission_Probabilities)}

    def set_circle_transition_probability_hmmer3(
            self,
            lTransition_Probabilities,
            final_state):
        """ Calculate transition probabilites between states that exist due to the
            circularity of the tandem repeat HMM.

        In HMMER3 models, no transition in circle from the last match state to the first
        match state are included.

        Values need to be found for the following transitions:
        m->m    m->i    m->d    i->m    i->i    d->m    d->d

        Compare ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf pp.84-85.

        Input: -ln(p) Output: Also -ln(p')

        When the HMM only has a single match state, i.e. lD == 1, the final match state
        equals the only match state, and no transition probabilities can be calculated
        via averaging over all other states. Therefore, the its transition probabilities
        remain unchanged.

        Args:
          lTransition_Probabilities (list of float): The first parameter.
          final_state (str): The second parameter.

        Returns:
          list of float: transition probabilities

        """

        if self.lD == 1:
            return lTransition_Probabilities

        # In self.hmmer, besides the probabilities for states 1, 2, 3, ... data for
        # "COMPO", 'letters', and 'id' are saved and should be ignored in the following.
        # Also, we ignore the transition probabilities for final_state, as our goal here
        # is to recalculate these based on the averages of the transition probabilities
        # from all other states.
        skip_states = ['COMPO', 'letters', final_state, 'id']

        # M->M, M->I, M->D: Use an average from all other transitions from
        # match states.
        log.debug("self.hmmer.keys(): %s", self.hmmer.items())
        log.debug("self.hmmer.items(): %s", self.hmmer.items())
        mean_M_M = np.mean([np.exp(-iP["transition"][0])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        mean_M_I = np.mean([np.exp(-iP["transition"][1])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        mean_M_D = np.mean([np.exp(-iP["transition"][2])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])

        sum_p = np.sum([mean_M_M, mean_M_I, mean_M_D])
        lTransition_Probabilities[
            :3] = [-np.log(i / sum_p) for i in [mean_M_M, mean_M_I, mean_M_D]]

        # I->M and I->I: Use an average from all other transitions from
        # insertion states.
        mean_I_M = np.mean([np.exp(-iP["transition"][3])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        mean_I_I = np.mean([np.exp(-iP["transition"][4])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        sum_p = np.sum([mean_I_M, mean_I_I])
        lTransition_Probabilities[
            3:5] = [-np.log(i / sum_p) for i in [mean_I_M, mean_I_I]]

        # D->M and D->D: Use an average from all other transitions from
        # deletion states.
        mean_D_M = np.mean([np.exp(-iP["transition"][5])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        mean_D_D = np.mean([np.exp(-iP["transition"][6])
                            for iState, iP in self.hmmer.items() if iState not in skip_states])
        sum_p = np.sum([mean_D_M, mean_D_D])
        lTransition_Probabilities[
            5:] = [-np.log(i / sum_p) for i in [mean_D_M, mean_D_D]]

        return lTransition_Probabilities

    def set_transition_probability_hmmer3(
            self,
            state_index,
            lTransition_Probabilities):
        """ Convert and assign transition probabilities from a HMMER3 file to HMM attribute `p_t`

        Set ``p_t`` of states to ``lTransition_Probabilities`` for the
        ``state_index`` th state

        In HMMER3 data, transition probabilities are -ln(p).
        Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)

        Args:
            state_index (int): E.g. 4
            lTransition_Probabilities (list of float): A list of transition probabilities
                as found in a single row for a certain state in HMMER3 files. E.g.
                ['0.00021', '8.85487', '9.57722', '0.61958', '0.77255', '0.48576', '0.95510']

        """

        lTransition_Probabilities = [-(i) * np.log10(np.exp(1))
                                     for i in lTransition_Probabilities]

        # Match state state_index -> Match state state_index + 1
        iM = self.match_states[state_index]
        iM_1 = self.match_states[(state_index + 1) % self.lD]
        self.p_t[iM][iM_1] = lTransition_Probabilities[0]

        # Match state state_index -> Insertion state state_index
        self.p_t[iM][
            self.insertion_states[state_index]] = lTransition_Probabilities[1]

        # Match state state_index -> Deletion state state_index + 1
        iD_1 = self.deletion_states[(state_index + 1) % self.lD]
        self.p_t[iM][iD_1] = lTransition_Probabilities[2]

        # Insertion state state_index ->  Match state state_index + 1
        iI = self.insertion_states[state_index]
        self.p_t[iI][iM_1] = lTransition_Probabilities[3]

        # Insertion state state_index ->  Insertion state state_index
        iI = self.insertion_states[state_index]
        self.p_t[iI][iI] = lTransition_Probabilities[4]

        # Deletion state state_index -> Match state state_index + 1
        iD = self.deletion_states[state_index]
        self.p_t[iD][iM_1] = lTransition_Probabilities[5]

        # Deletion state state_index -> Deletion state state_index + 1
        self.p_t[iD][iD_1] = lTransition_Probabilities[6]

    def get_direct_transition_probabilities_for_deletions(
            self,
            state_index,
            state):
        """ Calculate updated transition probabilities for HMMs without deletion states.

        Deletion states are non-emitting states. Therefore, they can be merged into the
        remaining emitting states. E.g. M -> D -> I would be equivalent to M -> I.
        The total transition probability from M to I would then be the product of the
        transition probabilities from M to D and from D to I (equivalent to the sum of
        the log10 probabilities.).

        This method translates all deletions following `state` with `state_index` into
        direct match state to match state transitions, and returns the new transition
        probabilities.

        Args:
            state_index (int): The first parameter.
            state (str): The second parameter.

        Returns:
            {state1: float}: The transition probability from `state` with `state_index` to
            `state1` is calculated for all (?) possible `state1`.
        """

        transition = {}
        p_Deletion_State = self.p_t[state][
            self.deletion_states[
                (state_index + 1) %
                self.lD]]
        for iState_index in range(1, self.lD):
            iDeleted_state_index = (state_index + iState_index) % self.lD
            iDeleted_state = self.deletion_states[iDeleted_state_index]
            iGoal_state_index = (state_index + iState_index + 1) % self.lD
            iGoal_state = self.match_states[iGoal_state_index]
            transition[iGoal_state] = p_Deletion_State + \
                self.p_t[iDeleted_state][iGoal_state]
            p_Deletion_State += self.p_t[iDeleted_state][
                self.deletion_states[iGoal_state_index]]

        return transition


############################### External parameters ######################

def hmmer3_emission_probabilities(hmmer_probabilities, letters, lMatch):
    '''
    Get emission probabilities from hmmer3 hmm file.
    In hmm file, emission probabilities are -ln(p).
    Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)

    Parameters (e.g.)::
        letters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
        lMatch = ['M'+str(i) for i in range(24)]

    Return format (pseudo code):
    [{iA: np.log10(p(iA,iM)) for iA in alphabet.keys()} for iM in lMatch]

    '''

    # Test: Is the number of match states in both models equal?
    if not len(hmmer_probabilities.keys()) - 2 == len(lMatch):
        print('Match states HMMER: {0} Match states local: {1}'.format(
            len(hmmer_probabilities.keys()) - 2, len(lMatch)))
        raise ValueError(
            'The number of match states in HMMer model and local model does not match')

    # Test: Are all `letters` represented in the HMMER HMM?
    if any(iL not in hmmer_probabilities['letters'] for iL in letters):
        missing = [
            iL for iL in letters if iL not in hmmer_probabilities['letters']]
        print('Missing representation in Hmmer File: {0}'.format(missing))
        raise ValueError(
            'Some letters in the local HMM are not represented in the HMMER HMM.')

    return [{iL: -(iP) * np.log10(np.exp(1)) for iL,
             iP in zip(hmmer_probabilities['letters'],
                       data['emissions']) if iL in letters} for key,
            data in hmmer_probabilities.items() if key not in ['letters',
                                                               'COMPO']]
