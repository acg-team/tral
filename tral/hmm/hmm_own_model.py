# (C) 2012-2015 Elke Schaper
# coding: utf-8

import math
import numpy as np
import scipy as sp
import scipy.stats
import scipy.special
import scipy.linalg
from scipy.optimize import fsolve, minimize_scalar
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.ERROR)

from tral.repeat import repeat
from tral.repeat.repeat_score import load_model
from tral.hmm import hmm_io

################################### HMM class ############################


class HMM:

    """ A cyclic HMM applicable to describe sequence Tandem Repeats """

    def __init__(
            self,
            tandem_repeat=None,
            prior_divergence=None,
            prior_indel_insertion=None,
            parameters=None):

        if tandem_repeat:
            if isinstance(tandem_repeat, str):
                self.HMM_from_file(
                    tandem_repeat,
                    accession=parameters,
                    prior_indel_insertion=prior_indel_insertion)
            else:
                self.HMM_from_TR(
                    tandem_repeat,
                    prior_divergence=prior_divergence,
                    prior_indel_insertion=prior_indel_insertion,
                    emission_file=parameters)
        else:
            self.HMM_example()

    def HMM_example(self):
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

    def HMM_from_TR_One_step(self, tandem_repeat, t):
        """ Build a HMM from a TR alignment given the ML divergence t.
            One step approach: Calc the posterior of the ancestral sequence and use it
            as the likelihood of the homologous repeat unit sequence. """
        # not implemented yet. For implementation, Use the first half of
        # self.HMM_from_TR()

    def HMM_from_file(self, hmm_file, accession, prior_indel_insertion):
        """ Load a HMM from hmm_file
            Store probabilities as log10arithms.

            Parameters:
            <hmm_file> = 'path//to/file.hmm'
            <accession> = 'PF00560'
            <prior_indel_insertion> = {'mu': 0.5, 'sigma_squared': 0.81}


        """

        bInsertion_state_emission_p_from_file = False
        bTerminal_state_emission_p_from_file = False

        hmmer_probabilities = hmm_io.read_HMMER(hmm_file, id=accession)
        dTranslate_States = {
            i: str(
                int(i) -
                1) for i in hmmer_probabilities.keys() if i not in [
                'COMPO',
                'letters']}
        hmmer_probabilities = {(i if i in ['COMPO', 'letters'] else dTranslate_States[
                                i]): j for i, j in hmmer_probabilities.items()}
        # WARNING: The order of keys in hmmer_probabilities might not be what
        # you think. You need to adapt your code at several instances (wherever
        # you assume the states to be ordered, e.g. zip(repeat_states,
        # hmmer_probabilities.keys))

        if not hmmer_probabilities:
            print("<hmmer_probabilities> were not found :(")
            return None

        self.hmmer = hmmer_probabilities

        # Warning! The following line is only appropriate, if there are only
        # "letters", "COMPO", and the match state keys in <hmmer_probabilities>
        l_effective = len(hmmer_probabilities.keys()) - 2
        # Assume: sequence_type = 'AA'
        Q, null_model_emission_p, alphabet = load_model('lg')

        # Initialise all HMM states to default value
        repeat_states, insert_states = self.initialise_HMM_structure(l_effective)

        # Store indel probabilities
        self.set_indel_probabilities(prior_indel_insertion, l_effective, default=True)

        # Store null model insertion probabilities for terminal and insertion
        # states
        if bTerminal_state_emission_p_from_file:
            null_model_emission_p = hmmer3_null_model_emission_probabilites(
                hmmer_probabilities=hmmer_probabilities)

        self.hmmer_probabilities = hmmer_probabilities

        # Store match state and insertion state emission probabilities
        match_state_emission_p = hmmer3_emission_probabilities(
            hmmer_probabilities=hmmer_probabilities,
            letters=list(
                alphabet.keys()),
            lMatch=repeat_states)
        if bInsertion_state_emission_p_from_file:
            insertion_state_emission_p = hmmer3_insertion_probabilities(
                hmmer_probabilities=hmmer_probabilities,
                letters=list(
                    alphabet.keys()),
                lMatch=insertion_states)
        else:
            insertion_state_emission_p = None
        self.set_emission_probabilities(
            alphabet,
            null_model_emission_p,
            match_state_emission_p,
            repeat_states,
            insertion_state_emission_p,
            insert_states)

    def HMM_from_TR(
            self,
            tandem_repeat,
            prior_divergence=None,
            prior_indel_insertion=None,
            emission_file=None):
        """ Build a HMM from a TR <tandem_repeat>.
            Two step approach: First, calc the posterior of the ancestral sequence
            Then, calculate the likelihood of the homologous repeat unit sequence.

            Parameters:
            <prior_divergence> is {'type': 'alpha'/'fixed_value', 'value': 4}
            <prior_indel_insertion> is {'mu': 0.5, 'sigma_squared': 0.81}

            Store probabilities as log10arithms."""

        # This will fail at current, as the divergence attribute is a dict.
        if prior_divergence is None:
            divergence = divergence_from_FP_simulations(tandem_repeat.l_effective)
        elif prior_divergence['type'] == 'alpha':
            divergence = divergence_from_FP_simulations(
                tandem_repeat.l_effective,
                alpha=prior_divergence['value'])
            # print(divergence)
            # print(prior_divergence['value'])
        else:
            divergence = prior_divergence['value']

        if prior_indel_insertion is not None and prior_indel_insertion[
                'type'] == 'adaptive':
            # Adapt the indel probabilities towards the divergences
            prior_indel_insertion['mu'] = divergence / \
                prior_indel_insertion['factor']
            prior_indel_insertion[
                'sigma_squared'] = prior_indel_insertion['mu'] / 5

        sequence_type = tandem_repeat.sequence_type
        msaTD = tandem_repeat.msaTD
        l_effective = tandem_repeat.l_effective
        n = tandem_repeat.n  # Test: Do you really need this variable?

        repeat_states, insert_states = self.initialise_HMM_structure(l_effective)

        # Transitions incorporating deletion or insertion states

        # We assume that deletions are not longer than the repeat unit, which would hence not be visible anymore.
        # We also assume that you cannot jump from "I" directly to "C", as is
        # would have been more likely to directly assume "C" from the foregoing
        # state.

        # # Apply the modulus in order to determine the true index of any state,
        # as you might have crossed the current tandem repeat unit's border:
        # index(i,l_effective) = i%l_effective

        insertion_lengths, deletion_lengths = tandem_repeat.gap_structure_HMM()
        self.set_indel_probabilities(
            prior_indel_insertion,
            l_effective,
            n=n,
            default=False,
            insertion_lengths=insertion_lengths,
            deletion_lengths=deletion_lengths)

        # Load model of sequence evolution parameters
        if sequence_type == 'AA':
            Q, eqFreq, alphabet = load_model('lg')
        else:
            Q, eqFreq, alphabet = load_model('tn93')

        # Calculate ML-emission probabilities for all match states
        # YOU might want to transfer Q,eqFreq,alphabet to this function???
        if not emission_file:
            likelihood_offspring = calculate_log10_offspring_likelihood(
                tandem_repeat,
                divergence)
        else:
            hmmer_probabilities = hmm_io.read_HMMER(emission_file)[0]
            likelihood_offspring = hmmer3_emission_probabilities(
                hmmer_probabilities=hmmer_probabilities,
                letters=list(
                    alphabet.keys()),
                lMatch=repeat_states)

        # Anpassen!
        self.set_emission_probabilities(
            alphabet=alphabet,
            null_model_emission_p=eqFreq,
            match_state_emission_p=likelihood_offspring,
            repeat_states=repeat_states)

    def initialise_HMM_structure(self, l_effective):
        '''
        Initialise all states
        Set transition probabilities to None/0
        Set initial probability for all states equal.
        '''

        # Build a HMM from these likelihoods.
        repeat_states = ["M{0}".format(str(i)) for i in range(l_effective)]
        insert_states = ["I{0}".format(str(i)) for i in range(l_effective)]
        self.states = ["N"] + repeat_states + insert_states + ["C"]
        logger.debug("HMM states: {0}".format(self.states))

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
        # - Transitions from "M_k" to "M_k+1" (for k<n) and "M_k" to "M_0" (for k==n)
        # - Transitions from "N" to "N"
        # - Transitions from "C" to "C"
        for i in range(l_effective):
            self.p_t["N"]["M{0}".format(index(i, l_effective))] = 0
            self.p_t["M{0}".format(index(i, l_effective))]["C"] = 0
            self.p_t[
                "M{0}".format(
                    index(
                        i,
                        l_effective))][
                "M{0}".format(
                    index(
                        i +
                        1,
                        l_effective))] = 0
        self.p_t["N"]["N"] = 0
        self.p_t["C"]["C"] = 0

        return repeat_states, insert_states

    def set_emission_probabilities(
            self,
            alphabet,
            null_model_emission_p,
            match_state_emission_p,
            repeat_states,
            insertion_emmission_p=None,
            insertion_states=None):

        # emissions
        self.emissions = alphabet.keys()
        logger.debug("Emissions: {0}".format(self.emissions))

        # emission probabilities

        # Set all emission probabilities to null model emission probabilities
        self.p_e = {
            iS: {
                letter: np.log10(
                    null_model_emission_p[number]) for letter,
                number in alphabet.items()} for iS in self.states}

        # Set match state probabilities to ML-emission probabilities
        for iRS, iEmission_p in zip(repeat_states, match_state_emission_p):
            self.p_e[iRS] = iEmission_p

        if insertion_emmission_p and insertion_states:
            print("Setting insertion state probabilities")
            for iIS, iEmission_p in zip(
                    insertion_states, insertion_emmission_p):
                self.p_e[iIS] = iEmission_p

        #logger.debug("Emission probabilities: {0}".format(self.p_e['M0']))

    def set_indel_probabilities(
            self,
            prior_indel_insertion,
            l_effective,
            n=2,
            default=True,
            insertion_lengths=None,
            deletion_lengths=None):

        if default:
            p_deletion_formation = calculate_log10_indel_probability(
                0,
                n,
                prior=prior_indel_insertion)
            p_insertion_formation = calculate_log10_indel_probability(
                0,
                n,
                prior=prior_indel_insertion)
            p_deletion_lengths = calculate_log10_probability_indel_lengths(
                [],
                l_effective -
                1,
                type='zipf')
            p_insertion_continued, p_insertion_stopped = calculate_log10_probability_indel_lengths(
                [], None, type='exponential')

        # First, calculate the transition probability of all states that only
        # depends on the gap structure of one site:
        for i in range(l_effective):

            if not default:
                p_deletion_formation = calculate_log10_indel_probability(
                    len(deletion_lengths[i]), n, prior=prior_indel_insertion)
                p_insertion_formation = calculate_log10_indel_probability(
                    len(insertion_lengths[i]), n, prior=prior_indel_insertion)
                p_deletion_lengths = calculate_log10_probability_indel_lengths(
                    deletion_lengths[i],
                    l_effective -
                    1,
                    type='zipf')
                p_insertion_continued, p_insertion_stopped = calculate_log10_probability_indel_lengths(
                    insertion_lengths[i], None, type='exponential')

            # In the current setup, moving from "N" to either an insertion, or a deletion state
            # is less likely than staying in "N". Thus, these probabilities are left at zero.
            #self.p_t["N"]["I{0}".format(str(i))] = p_insertion

            # Mi-1 -> Mi+l_effective (deletion(l_effective))
            for iC, p_deletion_length in enumerate(p_deletion_lengths):
                self.p_t[
                    "M{0}".format(
                        index(
                            i - 1,
                            l_effective))][
                    "M{0}".format(
                        index(
                            i + iC + 1,
                            l_effective))] = p_deletion_formation + p_deletion_length
                # Ii -> Mi + l_effective (insertionL) + deletion(l_effective)
                self.p_t[
                    "I{0}".format(
                        index(
                            i,
                            l_effective))][
                    "M{0}".format(
                        index(
                            i + iC + 1,
                            l_effective))] = p_insertion_stopped + p_deletion_formation + p_deletion_length

            # Mi-1 -> Ii (insertion)
            self.p_t[
                "M{0}".format(
                    index(
                        i - 1,
                        l_effective))][
                "I{0}".format(
                    index(
                        i,
                        l_effective))] = p_insertion_formation

            # Ii -> Ii (1-insertionL)
            self.p_t[
                "I{0}".format(
                    index(
                        i,
                        l_effective))][
                "I{0}".format(
                    index(
                        i,
                        l_effective))] = p_insertion_continued

            # Ii -> Mi (insertionL)
            self.p_t[
                "I{0}".format(
                    index(
                        i,
                        l_effective))][
                "M{0}".format(
                    index(
                        i,
                        l_effective))] = p_insertion_stopped

        # Second, calculate the transition probability of all states that depends on the gap structure of several sites.
        # These calcs use the information already gathered in self.p_t.
        # These are deletions, that are directly followed by insertions
        # (Case 1: deletion from match state. Case 2: deletion from insertion state.)
        for i in range(l_effective):
            for iC in range(1, l_effective):
                # Mi -> Ii+iC (deletion(iC-1) + insertion(i+iC))
                self.p_t[
                    "M{0}".format(
                        index(
                            i,
                            l_effective))][
                    "I{0}".format(
                        index(
                            i + iC,
                            l_effective))] = self.p_t[
                    "M{0}".format(
                        index(
                            i,
                            l_effective))][
                    "M{0}".format(
                        index(
                            i + iC,
                            l_effective))] + self.p_t[
                    "M{0}".format(
                        index(
                            i + iC - 1,
                            l_effective))][
                    "I{0}".format(
                        index(
                            i + iC,
                            l_effective))]

                # Ii -> Ii+l_effective insertion_stop(i) + deletion(i;l_effective) +
                # insertion_formation(i+l_effective+1)
                self.p_t["I{0}".format(index(i,
                                             l_effective))]["I{0}".format(index(i + iC,
                                                                       l_effective))] = self.p_t["I{0}".format(index(i,
                                                                                                            l_effective))]["M{0}".format(index(i,
                                                                                                                                      l_effective))] + self.p_t["M{0}".format(index(i - 1,
                                                                                                                                                                           l_effective))]["M{0}".format(index(i + iC,
                                                                                                                                                                                                     l_effective))] + self.p_t["M{0}".format(index(i + iC - 1,
                                                                                                                                                                                                                                          l_effective))]["I{0}".format(index(i + iC,
                                                                                                                                                                                                                                                                    l_effective))]

        logger.debug("Transition probabilities: {0}".format(self.p_t))


################################### Local TR class #######################
class TR:

    """ Fake interim TR class """

    def __init__(self):

        self.divergence = 0.6
        self.msaTD = ["AAA", "CCC", "GTG"]
        self.sequence_type = 'DNA'
        self.n = 3
        self.l_effective = 3

################################### MAP parameter estimation 3############


def divergence_from_FP_simulations(l, alpha=0.1):
    ''' Which HMM divergence sets the average FP-rate
        (i.e. the number of falsely positively assigend amino acid on either flanking side to a perfect TR)
        to in total alpha*l'''

    # alpha: In average, the number of falsely assigned AAs on both sides of the TR with should not extend l * alpha
    # Derivation:
    # FP-rate (average number of AAs predicted falsely as part of a perfect
    # TR) on one side of the TR as a function of the divergence.
    f = [
        0.50150,
        0.58600,
        0.70900,
        0.75125,
        0.89275,
        1.00900,
        1.25100,
        1.24150,
        1.56275,
        1.69600,
        1.88475,
        2.20025,
        2.33375,
        2.53500,
        2.99300,
        3.25900,
        3.77400,
        4.03300,
        4.26675,
        4.71800,
        5.30925,
        5.83325,
        6.31575,
        6.84000,
        7.93425,
        9.01900,
        9.72000,
        10.80925,
        12.26650,
        13.50750,
        14.29825,
        15.74225,
        18.66050,
        20.02575,
        21.54400,
        23.97300,
        26.93625,
        30.39400,
        32.78550,
        35.80375,
        38.98025,
        45.79925,
        48.48200,
        51.29100,
        56.74300,
        60.45925,
        64.44675,
        69.41900,
        80.02450,
        82.37550,
        88.47350,
        95.13425,
        98.00750,
        105.54975,
        108.85625,
        117.76750,
        125.69250,
        125.73725,
        137.40675,
        143.08225,
        148.28725,
        154.51425,
        161.23075,
        160.33875,
        172.00475,
        175.41525,
        180.44425,
        179.47350,
        190.14225,
        190.93625,
        194.03725,
        194.87475,
        203.35900,
        209.94350,
        210.63750,
        210.63950,
        216.69850,
        223.15025,
        219.22775,
        227.37700,
        227.09600,
        238.11900,
        235.79275,
        238.05250,
        242.59875,
        244.40425]
    # Sequence divergence associated with f.
    d = [(i / 10) + 0.5 for i in range(86)]
    # List <da> contains the linearly weighted averaged divergences of the two
    # values of the FP-rate that enclose l*alpha:
    da = []
    before = (f[0], d[0])
    for iF, iD in zip(f[1:], d[1:]):
        f_ok = (alpha * (len(da) + 1)) / 2
        while iF > f_ok:
            da.append(
                ((iF - f_ok) * before[1] + (f_ok - before[0]) * iD) / (iF - before[0]))
            f_ok = (alpha * (len(da) + 1)) / 2
        before = (iF, iD)

    if l > len(da):
        return 6
    else:
        return da[l - 1]


def index(i, iMax):
    return str(i % iMax)


def calculate_log10_offspring_likelihood(tandem_repeat, divergence=None):

    sequence_type = tandem_repeat.sequence_type
    msaTD = tandem_repeat.msaTD
    l_effective = tandem_repeat.l_effective
    n = tandem_repeat.n

    # Calculate posterior probabilities of ancestral states from the TR
    # alignment

    # Load model of sequence evolution parameters
    if sequence_type == 'AA':
        Q, eqFreq, alphabet = load_model('lg')
    else:
        Q, eqFreq, alphabet = load_model('tn93')
    P = scipy.linalg.expm(Q * divergence)
    alphabet_reverse = {i: j for j, i in alphabet.items()}

    # Initialise the count matrix (if not already existent)
    try:
        tandem_repeat.msaTDN
    except:
        msaTDN_temp = [[alphabet[symbol] for symbol in column if symbol not in [
            '-', 'X', '*', 'U']] for column in tandem_repeat.msaTD]
        tandem_repeat.msaTDN = []
        for column in msaTDN_temp:
            my_diag = {i: 0 for i in range(len(alphabet.keys()))}
            for iN in column:
                my_diag[iN] = my_diag[iN] + 1
            tandem_repeat.msaTDN.append(
                np.array([j for j in my_diag.values()]))

    # Initialise results matrix
    posterior_ancestor = []

    # Calculate the posterior probabilities per column
    for column in tandem_repeat.msaTDN:
        posterior = {
            i: iJ *
            np.product(
                np.power(
                    P[i],
                    column)) for i,
            iJ in enumerate(eqFreq)}
        normalising_factor = sum(posterior.values())
        for key, value in posterior.items():
            posterior[key] = value / normalising_factor
        posterior_ancestor.append(posterior)

    # Calculate the likelihood of evolved TR units with the given divergence
    # from these ancestral states
    likelihood_offspring = [
        {
            iA: np.log10(
                sum(
                    posterior[iAncestor] *
                    P[iAncestor][i] for iAncestor in alphabet.values())) for iA,
            i in alphabet.items()} for posterior in posterior_ancestor]
    return likelihood_offspring


def loglikelihood_substitution(t, Q, eqFreq, alphabet, tandem_repeat):
    ''' calculate the likelihood of a repeat assuming a star tree and a sequence model
        defined by Q, t, eqFreq and alphabet. '''
    P = scipy.linalg.expm(Q * t)

    # msaTDN is a list of tandem_repeat.n numpy arrays.
    # Each numpy array represents the elements of the alphabet. The entries
    # are the count of an alphabet element in the respective tandem_repeat
    # column
    try:
        tandem_repeat.msaTDN
    except:
        msaTDN_temp = [[alphabet[symbol] for symbol in column if symbol not in [
            '-', 'X', '*', 'U']] for column in tandem_repeat.msaTD]
        tandem_repeat.msaTDN = []
        for column in msaTDN_temp:
            my_diag = {i: 0 for i in range(len(alphabet.keys()))}
            for iN in column:
                my_diag[iN] = my_diag[iN] + 1
            tandem_repeat.msaTDN.append(
                np.array([j for j in my_diag.values()]))

    # First sum: Over all sites
    # Second sum: Over all possible ancestral characters
    # Third product: Over all repeat units
    likelihood = sum(
        np.log(
            np.sum(
                iJ *
                np.product(
                    np.power(
                        P[i],
                        column)) for i,
                iJ in enumerate(eqFreq))) for column in tandem_repeat.msaTDN)
    return likelihood


def calculate_log10_probability_indel_lengths(
        indel_lengths,
        indel_length_max,
        type='zipf',
        prior=None):
    ''' calculate the probability of indels of different length until a maximum length <indel_length_max>,
    assuming that either

    type == 'zipf'
    indel lengths are distributed following a Zipfian distribution with parameter <indel_zipf>.
    (Compare Fletcher,W. and Yang,Z. (2009) INDELible: a flexible simulator of biological sequence evolution. Mol Biol Evol, 26, 1879–1888.
    In their publication, <indel_zipf> is denoted as <a>, the gap lengths as <u>.)

    or
    type == 'exponential'
    indels lengths are distributed following a (1-alpha)*(alpha^(l-1)) distribution.
    [Is this model already used? Probably yes. How can i cite it?]

    A good choice for <indel_length_max> when building HMMs for tandem repeats is <l_effective-1>, as a deletion, that
    includes a complete repeat unit, will not be visible in the tandem repeat unit alignment, and could therefore
    be ignored.

    Return the probabilities as a list.
    '''

    if 'zipf' == type:
        indel_zipf = calculate_MAP_Indel_length_Zipfian_factor(
            indel_lengths,
            prior)
        zeta_factor = np.log10(1 / scipy.special.zeta(indel_zipf, 1))
        return [
            np.log10(
                ((iL + 1) ** -indel_zipf)) + zeta_factor for iL in range(indel_length_max)]

    elif 'exponential' == type:
        alpha = calculate_MAP_Indel_length_exponential_factor(
            indel_lengths,
            prior)
        return np.log10(alpha), np.log10(1 - alpha)

    else:
        logger.error(
            "I am not aware of indel length model type: {0}.".format(type))


def calculate_MAP_Indel_length_exponential_factor(indel_lengths, prior=None):
    ''' calculate the MAP Exponential decay constant <alpha> that determines the distribution of
    indel lengths. Assume a Gaussian prior. Input is a list of <indel_lengths> for a particular
    column.

    The MAP is calculated numerically, as no nice analytical solution has been found so far [power function of grade 3 needs to be solved:
    0 == prior['sigma_squared']*(sum(indel_lengths) - len(indel_lengths)) - alpha*(prior['sigma_squared']*sum(indel_lengths) + prior['mu']) + (alpha**2)*(1+prior['mu']) - alpha**3
    '''

    if prior is None:
        prior = {'mu': 0.8, 'sigma_squared': 0.1}

    if len(indel_lengths) == 0:
        return prior['mu']

    else:
        # As we are only interested in the maximum, we can leave out factors (e.g. * 1/(math.sqrt(2*math.pi * prior['sigma_squared'])) for the Gaussian)
        # There is an analytic solution. Compare "Gap it!" or "gaps.nb"
        posterior = lambda alpha: - np.prod([(1 - alpha) * (alpha ** (lIndel - 1))
                                             for lIndel in indel_lengths]) * math.exp(-0.5 * ((alpha - prior['mu']) ** 2) / prior['sigma_squared'])
        res = minimize_scalar(posterior, bounds=(0, 1), method='bounded')
        return res.x


def calculate_MAP_Indel_length_Zipfian_factor(indel_lengths, prior=None):
    ''' calculate the MAP Zipfian constant <indel_zipf> that determines the distribution of
    indel lengths. Assume a Gaussian prior. Input is a list of <indel_lengths> for a particular
    column.
    The probability distribution of indel lengths is assumed to follow the Zipfian distribution
    (Compare Fletcher,W. and Yang,Z. (2009) INDELible: a flexible simulator of biological sequence evolution. Mol Biol Evol, 26, 1879–1888.
    In their publication, <indel_zipf> is denoted as <a>, the gap lengths as <u>.)

    The MAP is calculated numerically, as there did not seem to be a nice analytical solution.
    '''

    if prior is None:
        prior = {'mu': 1.821, 'sigma_squared': 1}

    if len(indel_lengths) == 0:
        return prior['mu']

    else:
        # The posterior is the likelihood * prior * scaling.
        # As we are only interested in the maximum, we leave out all factors in
        # the following equation for the posterior:
        posterior = lambda indel_zipf: - (scipy.special.zeta(indel_zipf, 1) ** - len(indel_lengths)) * np.prod(
            [lIndel ** -indel_zipf for lIndel in indel_lengths]) * math.exp(-0.5 * ((indel_zipf - prior['mu']) ** 2) / prior['sigma_squared'])
        res = minimize_scalar(posterior, method='brent')
        return res.x


def calculate_log10_indel_probability(nIndels, n, prior=None):
    ''' calculate the probabilty of an indel per site and per repeat unit
    from the MAP estimate of the indel rate, given
    - the <divergence> of the repeat units
    - a <prior> on the indel rates (normal distributed)
    - the number of indels <nIndels> in this column
    - the total length <n> of this column
    '''

    indel_rate_MAP = calculate_MAP_indel_rate(nIndels, n, prior)
    return np.log10(1 - math.exp(-indel_rate_MAP))


def calculate_MAP_indel_rate(nIndels, n, prior=None):
    ''' calculate the MAP indel_rate per repeat unit and site given
    - an assumed prior distribution of the indel rate described by the dict <prior>
    - the likelihood of the the observed number of indels <nIndels>.

    posterior ~ likelihood * prior
    Derive the posterior with respect to the indel_rate
    Set the result equal to 0
    Solve with respect to the MAP of the indel_rate.

    For nIndels in {0,1} it is possible to derive the MAP indel_rate analytically.
    For all other cases, the maximum of the posterior distribution needs to be found in good approximation algorithmically.

    ln(posterior) ~ ln(likelihood) + ln(prior)
    '''

    if prior is None:
        prior = {'mu': 0.01, 'sigma_squared': 0.001}

    if nIndels == 0:
        # This max is chosen absolutely randomly. do you have any better idea?
        # (Before: constant 0.001)
        return max(prior['mu'] / 100, prior['mu'] - prior['sigma_squared'] * n)
    else:
        # Calculate the roots numerically
        res = fsolve(
            derivative_log_posterior,
            prior['mu'],
            args=(
                prior,
                n,
                nIndels))
        return res[0]


def derivative_log_posterior(indel_rate, prior, n, nIndels):
    ''' assuming that on each site on each repeat unit, an in/del event has occured, or not.
    Thus, the likelihood of a certain number of in/del events follows the binomial distribution.
    The prior is assumed to be a gaussian described by prior.

    Return the derivative with respect to <indel_rate> of the natural logarithm of the posterior
    probability of <nIndels>, calculated as the product of likelihood and prior. The
    normalisation factor is ignored for the moment.'''
    return indel_rate + (prior['mu'] - indel_rate - n * prior['sigma_squared']) * \
        math.exp(-indel_rate) - prior['mu'] + (n - nIndels) * prior['sigma_squared']

############################### External parameters ######################


def hmmer3_emission_probabilities(hmmer_probabilities, letters, lMatch):
    '''
    Get emission probabilities from hmmer3 hmm file.
    In hmm file, emission probabilities are -ln(p).
    Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)

    Parameters (e.g.):
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

    # Test: Are all <letters> represented in the HMMER HMM?
    if any(iL not in hmmer_probabilities['letters'] for iL in letters):
        missing = [
            iL for iL in letters if iL not in hmmer_probabilities['letters']]
        print('Missing representation in Hmmer File: {0}'.format(missing))
        raise ValueError(
            'Some letters in the local HMM are not represented in the HMMER HMM.')

    return [
        {
            iL: -
            float(iP) *
            np.log10(
                np.exp(1)) for iL,
            iP in zip(
                hmmer_probabilities['letters'],
                hmmer_probabilities[
                    iMatch[
                        1:]]['emissions']) if iL in letters} for iMatch in lMatch]


##################################### Tests ##############################


def test():
    ''' To be implemented... '''
    #tandem_repeat = ...
    divergence = 0
    calculate_log10_offspring_likelihood(tandem_repeat, divergence)

##################################### Main ###############################


def main():
    my_TR = repeat_info.Repeat(
        begin=0,
        msa=[
            'A-G',
            'ACG',
            'ACG'],
        sequence_type='DNA')
    #my_TR = TR()
    my_HMM = HMM(my_TR, divergence=0.1)

    #print(calculate_indel_rate(1, 6, 1.5))

    # print(calculate_probability_indel_lengths(0.366869,6,'exponential'))

if __name__ == "__main__":
    main()
