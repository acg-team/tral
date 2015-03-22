# (C) 2012-2015 Elke Schaper
from collections import defaultdict
import logging
import numpy as np
import operator
import os
import re

from tral import configuration
from tral.paths import *
from tral.repeat import repeat

log = logging.getLogger(__name__)

c = configuration.Configuration.Instance()
config = c.config


def viterbi(hmm, emission):
    """ Calculate the most probable sequence of states given a sequence of emissions
    and a HMM using the Viterbi algorithm

    Get local copies of all variables.
    All probabilities must be given as logarithms
    Replace Selenocysteine (U) with Cysteine (C), replace Pyrrolysine (O) with Lysine (K)
    [Reason: Seldom AAs that are not part of standard amino acid substitution models.]

    Args:
        hmm (hmm): An instance of the HMM class.
        emission (sequence): An instance of the Sequence class.

    Returns:
        The most likely sequence of hmm states to emit the sequence in the
        form of a list of str.

    .. todo:: Adapt docstrings to refactored Viterbi -> Viterbi_path classes.
    .. todo:: Check: Do you need local copies of all variables?
    .. todo:: Do the functions related to viterbi need to be summarized (e.g. in one
         class?) How do they relate to the Sequence class, or the HMM class?
    """

    if len(emission) / \
            hmm.l_effective < float(config['filter']['basic']['dict']['n_effective']['threshold']):
        logging.info(
            "Skip the HMM as it is too long ({}) for this sequence ({}) according to the filter criterion min n_effective ({}).".format(
                hmm.l_effective,
                len(emission),
                config['filter']['basic']['dict']['n_effective']))
        return None
    if hmm.l_effective > float(config['hmm']['l_effective_max']):
        logging.info(
            "Skip the HMM as it is too long ({}) according to the filter criterion max hmm.l_effective ({}).".format(
                hmm.l_effective,
                config['hmm']['l_effective_max']))
        return None

    states = hmm.states
    p_0 = {iS: value for iS, value in hmm.p_0.items()}
    p_e = {
        iS: {
            iE: value for iE,
            value in emission.items()} for iS,
        emission in hmm.p_e.items()}
    p_t = {
        iS: {
            iT: value for iT,
            value in transition.items()} for iS,
        transition in hmm.p_t.items()}

    if any([(iS not in config['lAll_amino_acid']) for iS in emission]):
        raise Exception(
            "There is an unknown amino acid in:\n {}\n".format(emission))

    # In case there are ambiguous amino acids in the sequence, calculate the expected frequencies of the AAs that they could stand for
    # from the emission frequencies of the hmm in the neutral state "N".
    dAmbiguous_local = {}
    for iA in config['dAmbiguous_amino_acid'].keys():
        dAmbiguous_local[iA] = {}
        if iA in emission:
            total = np.log10(sum(
                10 ** p_e['N'][iAmbiguous] for iAmbiguous in config['dAmbiguous_amino_acid'][iA]))
            for iAmbiguous in config['dAmbiguous_amino_acid'][iA]:
                dAmbiguous_local[iA][iAmbiguous] = p_e['N'][iAmbiguous] - total

    logging.debug("p_0: {0}".format(p_0))
    logging.debug("p_e: {0}".format(p_e))
    logging.debug("p_t: {0}".format(p_t))

    # Initialisation of the probabilities on the first emitted character
    if emission[0] in dAmbiguous_local:
        path = {iS: {'probability': p_0[iS], 'path': [iS]} for iS in states}
        for iS in states:
            # Calculate the average emission probability of ambiguity chars.
            # The numerical trick to calculate the average log values in high
            # precision is reused in the next for-loop and described there.
            lP_emission = [
                dAmbiguous_local[
                    emission[0]][iAmbiguous] +
                p_e[iS][iAmbiguous] for iAmbiguous in dAmbiguous_local[
                    emission[0]]]
            max_p = max(lP_emission)
            path[iS][
                'probability'] += np.log10(sum([10 ** (i - max_p) for i in lP_emission])) + max_p

    else:
        path = {
            iS: {
                'probability': p_0[iS] +
                p_e[iS][
                    emission[0]],
                'path': [iS]} for iS in states}
    logging.debug("Path: {0}".format(path))

    # Viterbi on all remaining emitted characters
    for iE in emission[1:]:
        logging.debug("Emitted: {0}".format(iE))

        # Determine next most probable state and its probability ...
        for iS in states:
            p = {}
            for iFormer in states:
                if iE in dAmbiguous_local:
                    # Calculate the probability of being in state iS and emitting iAmbiguous for any of the AAs that the ambiguous iE stands for.
                    # Then, average over these probabilities (taking into
                    # account the background frequencies of all iAmbiguous)
                    dP_Former = {
                        iAmbiguous: probability_of_the_former_state(
                            iFormer,
                            iS,
                            iAmbiguous,
                            p_e,
                            p_t,
                            path,
                            p) for iAmbiguous in dAmbiguous_local[iE].keys()}
                    dP_Former = {
                        iAmbiguous: j for iAmbiguous,
                        j in dP_Former.items() if j}
                    if len(dP_Former) > 0:
                        # Next, we need to calculated a weighted average over log10 probabilities.
                        # For this purpose, we need to transform the log10ps back to ps.
                        # However, the ps might have too small values. We apply a numerical trick:
                        # Instead of sum(w(i)p(i)) = sum ( 10 ** ( log(w(i)) + log(p(i)) ) ) we calculate
                        # sum(w(i)p(i)) = sum ( 10 ** ( log(w(i)) + log(p(i)) + k ) ) / (10 ** k), which can easily be shown to be equivalent.
                        # For (-k), we choose the maximum value in log(w(i)) +
                        # log(p(i)), as it corresponds to the maximum and thus
                        # most pronounced probability of the sum.
                        lP_Former = [
                            i +
                            dAmbiguous_local[iE][iAmbiguous] for iAmbiguous,
                            i in dP_Former.items()]
                        max_p = max(lP_Former)
                        p[iFormer] = np.log10(
                            sum([10 ** (i - max_p) for i in lP_Former])) + max_p
                else:
                    p_Former = probability_of_the_former_state(
                        iFormer,
                        iS,
                        iE,
                        p_e,
                        p_t,
                        path,
                        p)
                    if p_Former:
                        p[iFormer] = p_Former

            # This error will occur when it is not possible to enter a state iS at emission iE. This can happen for more complex HMMs, but should not happen for the simple
            # circular HMMs that we consider here. It means that state iS
            # cannot be reached by emission iE. Either iS cannot emit iE, or iS
            # cannot be reached by transitions at this point.
            if p == {}:
                logging.error("The dict p is empty. This should never happen.")
            # Detect the most likely path.
            next_state, p_next_state = max(
                p.items(), key=operator.itemgetter(1))
            path[iS]['probability_next'] = p_next_state
            path[iS]['path_next'] = path[next_state]['path'][:]
            path[iS]['path_next'].append(iS)

        # ... and update the path accordingly
        for iV in path.values():
            iV['path'] = iV['path_next'][:]
            iV['probability'] = iV['probability_next']
        #logging.debug("Path: {0}".format(path))

    # Which overall path is the most likely?
    path_summary = {iS: path[iS]['probability'] for iS in states}
    most_likely_terminal_state, p_most_likely_path = max(
        path_summary.items(), key=operator.itemgetter(1))
    most_likely_path = path[most_likely_terminal_state]['path']
    logging.debug(
        "The most likely path is {0}. It has a score of {1}.".format(
            most_likely_path,
            p_most_likely_path))

    return(most_likely_path)


def probability_of_the_former_state(iFormer, iS, iE, p_e, p_t, path, p):
    ''' Calculate the probability of iFormer, given iS and iE together with the dicts of emission and transition probabilities.'''
    # Exclude paths with a zero probability (i.e. a log(probability) set to
    # None)
    if iS in p_t[iFormer] and None is not p_t[iFormer][iS] and iE in p_e[
            iS] and None is not p_e[iS][iE] and None is not path[iFormer]['probability']:
        return path[iFormer]['probability'] + p_t[iFormer][iS] + p_e[iS][iE]

    return None


def distance_index(i, j, length):
    """ Helper function to calculate the distance between two indices in a circular HMM.

    Args:
        i (int): first index
        j (int): second index
        length (int): length of the HMM

    Returns:
        int: The distance between the indices ``i`` and ``j``: As the HMM
        is circular ``j`` may have a smaller value than ``i``, even though
        ``j`` is ahead of ``i``.

    .. todo:: May be replaced by a simple mod, as distance_index is used only once in
        the code at current.
    """
    if j > i:
        return j - i
    else:
        return length + j - i


def hmm_path_to_maximal_complete_tandem_repeat_units(
        lSequence,
        lPath,
        l_effective,
        alpha=None):
    """ Convert several viterbi paths of a hmm on several sequences into the corresponding hmm units.

    Be ungreedy: Start from the last index in the cluster of all start state and end state indices.

    Only integrate repeat units that are at least alpha complete (be it before of after)
    If you prefer absolute number of characters to filter which repeat units are used, and which not,
    change this in two occurrences of alpha.

    Assume that all states are counted starting on 1.

    Args:
        lSequence (list of str): A list of sequences.
        lPath (list of list of str): A list of Viterbi paths
        l_effective (int): length of the HMM used to create the Viterbi paths.
        alpha (Float): alpha element [0,1].

    Returns:
        A list of multiple sequence alignments (MSAs) in the form of a list of
        list of str, e.g. ``[['ATAILC', 'ATAILC', 'ACALKG'], ...]``.

    Raises:
        Exception: If alpha is not in [0,1].

    .. todo:: Should the sequences be sequence instances instead of just strings?
    .. todo:: Check example for returns.

    """

    if not alpha:
        alpha = 0.6

    # Find the first match state index used by all paths (add None if None)
    lTerminal_Indices = []
    lTRPresent = []
    for path in lPath:
        for iP in path:
            if iP.startswith("M"):
                lTerminal_Indices.append(int(iP[1:]))
                lTRPresent.append(True)
                break
        else:
            lTRPresent.append(False)

        for iP in path[::-1]:
            if iP.startswith("M"):
                lTerminal_Indices.append((int(iP[1:])) % l_effective + 1)
                break

    # Find the max distance between two indices
    lUsed_indices = sorted(
        [i for i in set(lTerminal_Indices) if i is not None])

    if len(lUsed_indices) == 0:
        return None
    elif len(lUsed_indices) == 1:
        start_index = lUsed_indices[0]
    else:
        distances = [i - j for j,
                     i in zip(lUsed_indices[:-1],
                              lUsed_indices[1:])] + [l_effective - lUsed_indices[-1] + lUsed_indices[0]]
        max_distance_index, max_distance = max(
            enumerate(distances), key=operator.itemgetter(1))

        # Define the start_index
        start_index = lUsed_indices[
            (max_distance_index + 1) %
            len(lUsed_indices)]

    # Get all TR units according to shift.
    lMSA = []

    for ibPresent, iSeq, iPath in zip(lTRPresent, lSequence, lPath):

        if not ibPresent:
            lMSA.append(None)
            continue

        msa = []
        # Intialise the last used index and the first msa_unit
        current_index = 1
        msa_unit = ''
        for iS, iP in zip(iSeq, iPath):
            if iP == 'N':
                continue
            elif iP == 'C':
                break
            match_state_index = (int(iP[1:]) + (l_effective - start_index)) % l_effective + 1
            if match_state_index >= current_index:
                # We are staying within the same repeat unit.
                msa_unit += iS
            else:
                # We are starting a new repeat unit.
                msa.append(msa_unit)
                msa_unit = iS
            current_index = match_state_index

        msa.append(msa_unit)

        try:
            if (len(msa[-1]) < alpha * l_effective):
                msa = msa[:-1]
            if (len(msa[0]) < alpha * l_effective):
                msa = msa[1:]
            lMSA.append(msa)
        except:
            lMSA.append([])

    return lMSA


def hmm_path_to_non_aligned_tandem_repeat_units(sequence, path, l_effective):
    """ Convert a viterbi <path> of a hmm of length <l_effective> on <sequence> into the corresponding tandem repeat

    Extract the tandem repeat alignment from a sequence given a Viterbi path.
    Ignore the alignment information in the Viterbi path. For example, all emissions
    labelled with M1 (match state 1) align according to the HMM. However, this method does
    not use this information. Therefore, for example the first characters in the repeat
    units do not necessarily align, and the repeat units are not necessarily of same length.

    Assume that all states are counted starting on 1.

    Args:
        sequence (str): A sequence as a string.
        lPath (list of list of str): A list of Viterbi paths
        l_effective (int): length of the HMM used to create the Viterbi paths.

    Returns:
        A multiple sequence alignment (MSA) created from the most likely path
        along the hmm in the form of a list of str.

    .. todo:: Should the sequence be a sequence instance instead of just a string?

    """

    begin = path.count('N')
    # If no repeat was found, return None.
    if begin == len(sequence) or path.count('C') == len(sequence):
        return None

    splitter = re.compile("(\w)(\d+)")
    mapping = [
        ((splitter.match(iP).group(1), int(
            splitter.match(iP).group(2))), iS) for iS, iP in zip(
            sequence[
                begin:], path[
                    begin:]) if iP != 'C']

    shift = mapping[0][0][1]
    index_shift = ["Empty"] + [(i + l_effective + 1 - shift) % l_effective for i in range(l_effective)]
    logging.debug("The tandem repeat is shifted by: {0}".format(str(shift)))

    repeat_msa = []
    repeat_unit = mapping[0][1]
    last_used_index = 0

    if l_effective == 1:
        for iM, iS in mapping[1:]:
            if iM[0] == "I":
                repeat_unit += iS
            else:
                repeat_msa.append(repeat_unit)
                repeat_unit = iS

    else:
        for iM, iS in mapping[1:]:
            #print("{} {}".format(iM,iS))
            if index_shift[iM[1]] >= last_used_index:
                repeat_unit += iS
            else:
                repeat_msa.append(repeat_unit)
                repeat_unit = iS
            last_used_index = index_shift[iM[1]]

    repeat_msa.append(repeat_unit)

    return repeat_msa


def hmm_path_to_aligned_tandem_repeat_units(sequence, most_likely_path, l_effective,
                                            translate=False):
    """Convert a viterbi path in an hmm of length ``l_effective`` on the sequence into a
    corresponding tandem repeat.

    Extract the tandem repeat alignment from a sequence given a Viterbi path.
    Use alignment information in the Viterbi path. For example, all emissions
    labelled with M1 (match state 1) align according to the HMM. Insert gaps
    for insertions and deletions accordingly.
    Thus, for example the first characters in the repeat units do necessarily
    align, albeit some of them may be gaps. Also, all repeat units are
    necessarily of same length.

    Assume that all states are counted starting on 0, unless the ``translate``
    flag is set.

    Args:
        sequence (str): A sequence as a string.
        most_likely_path (list of str): a Viterbi path.
        l_effective (int): length of the HMM used to create the Viterbi paths.
        translate (bool): This function assumes that HMM states are enumerated
            starting on 0.
            If the HMM states are enumerated starting on 1, set this flag for
            transformation.

    Returns:
        The function returns a tuple consisting of 3 values.
        The tuple contains:

        * ``msa``: A repeat instance.
        * ``begin``: The start index of the repeat on the sequence.
        * ``shift``: The index of the HMM where the cut between repeat units is set.


    .. warning:: [None].
    .. todo:: Should the sequence be a sequence instance instead of just a string?
    .. todo:: Check: How is the returned `begin` defined? Starting counting on 0 or 1?
        Is it the index of the last flanking character, or the first repeat character?
    .. todo:: Can we update this function, e.g. to not assume that HMM states start on 0?
    .. todo:: Check the docstring, reformat returns.
    """

    begin = most_likely_path.count('N')
    # If no repeat was found, return None.
    if begin == len(sequence) or most_likely_path.count('C') == len(sequence):
        return None

    if translate:
        # In the input data, the states are counted starting on 1. Subtract 1.
        most_likely_path = [
            i if i in ["C", "N"] else i[0] + str(int(i[1:]) - 1) for i in most_likely_path]

    splitter = re.compile("(\w)(\d+)")

    mapping = [
        ((splitter.match(iP).group(1), int(
            splitter.match(iP).group(2))), iS) for iS, iP in zip(
            sequence[
                begin:], most_likely_path[
                    begin:]) if iP != 'C']

    shift = mapping[0][0][1]
    index_shift = [(i + l_effective - shift) % l_effective for i in range(l_effective)]
    logging.debug("The tandem repeat is shifted by: {0}".format(str(shift)))

    insertions = []
    repeat_text = mapping[0][1]
    insertions.append(defaultdict(str))
    max_used_index_M = 0
    max_used_index_I = -1
    last_used_index = shift

    for iM, iS in mapping[1:]:

        logging.debug(
            "Iteration iS: {0}, iM: {1}, last_used_index: {2}, max_used_index_M: {3}, max_used_index_I: {4}".format(
                str(iS),
                str(iM),
                str(last_used_index),
                str(max_used_index_M),
                str(max_used_index_I)))

        # If we have entered a new repeat unit, add a new element to <insertions>
        # (Including some index magic, e.g. the insertion state index is shifted by one (lowered))
        if iM[0] == "M" and (
            index_shift[
                iM[1]] <= max_used_index_M or index_shift[
                iM[1]] <= max_used_index_I):
            insertions.append(defaultdict(str))
            max_used_index_M = index_shift[iM[1]]
            max_used_index_I = index_shift[iM[1]] - 1

        elif iM[0] == "I" and (index_shift[iM[1] - 1] < max_used_index_I or index_shift[iM[1] - 1] < max_used_index_M):
            insertions.append(defaultdict(str))
            max_used_index_I = index_shift[iM[1] - 1]
            max_used_index_M = index_shift[iM[1] - 1]

        if iM[0] == "M":
            max_used_index_M = index_shift[iM[1]]
        else:
            max_used_index_I = index_shift[iM[1] - 1]

        # Save match state and deletion information
        if iM[0] == "M":
            n_deletions = distance_index(last_used_index, iM[1], l_effective) - 1
            repeat_text += "-" * n_deletions + iS
            last_used_index = iM[1]

        else:
            insertions[-1][index_shift[iM[1]]] += iS

    msa = [repeat_text[i:i + l_effective] for i in range(0, len(repeat_text), l_effective)]
    msa[-1] += "-" * (l_effective - len(msa[-1]))
    msaT = ["".join(c) for c in zip(*msa)]
    n = len(msa)
    logging.debug("This tandem repeat has {0} repeat units.".format(str(n)))
    logging.debug(
        "These insertions were detected: {0}.".format(
            str(insertions)))

    msaTemp = []
    # For each site ...
    for i in range(l_effective):
        msaTemp.extend(msaT[i:i + 1])
        # ... for each tandem repeat unit, check whether there are insertions.
        for j, row in enumerate(insertions):
            # Backtranslate the index used in insertions (shift by one to the
            # right (increase))
            if (i + 1) % l_effective in row:
                for iS in row[(i + 1) % l_effective]:
                    msaTemp.append("-" * j + iS + "-" * (n - j - 1))
    msa = ["".join(c) for c in zip(*msaTemp)]

    return msa, begin, shift
