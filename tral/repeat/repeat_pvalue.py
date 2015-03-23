# (C) 2015 Elke Schaper

import os
import bisect
from os.path import join
import csv
import logging
import numpy as np
import scipy as sp

from tral.paths import DATA_DIR

LOG = logging.getLogger(__name__)

path_score = join(DATA_DIR, 'pvalue')

# ######################### REPEAT SCORE P-VALUE CALCULATION FUNCTIONS ########
# #################### phylo & entropy ########################################


def empirical_list(l, n, sequence_type='AA', score_type='phylo'):
    """Load and return a numpy list with 10,000 empirical values of the user defined
        distribution from an external file.


    If no distribution is available for a given ``l`` and ``n``, the closest values for
    ``l`` and ``n`` are chosen instead.

    For the standard model score "phylo", currently values are available until

        * ``l_max`` = 99
        * ``n_max`` = 49
        * ``total_repeat_length_max`` = 1000

    For ``l`` <=  ``l_max`` and ``n`` <=  ``n_max`` and ``n*l`` <=  ``total_repeat_length_max``
    all values are available.

    For other values, we use the closest distribution, assuming that values change little
    (heuristic assumption!).

    Args:
        l (int): The length of the repeat unit.
        n (int): The number of repeat units.
        sequence_type (str): The type of the sequence: either "AA" or "DNA".
        score_type: The model of repeat evolution used. E.g. "phylo".

    Returns:
        A numpy list of length 10.000


    .. todo:: Define "phylo" model.
    .. todo:: Perhaps add distributions beyond l_max, n_max.
    """

    l_max = 99
    n_max = 49
    if score_type == 'entropy':
        l_max = 20
        n_max = 5

    total_repeat_length_max = 1000
    l = min(l, l_max)
    n = min(n, n_max)
    if l * n > total_repeat_length_max:
        LOG.debug(
            "l: %d and n: %d are bigger than the total_repeat_length_max: %d" %
            (l, n, total_repeat_length_max))
        # Dirty little hack: as we do not have data for this pair of l and n,
        # apply nearest data file.
        while l * n > total_repeat_length_max:
            if l > n:
                l -= 1
            else:
                n -= 1

    file = join(path_score, sequence_type, score_type,
                str(l) + '_' + str(n) + '.npz')
    if not os.path.isfile(file):
        raise ValueError("complete pdf file %s does not exist!" % file)

    myEmpiricalList = np.load(file)
    # It is absolutely necessary to close the numpy filehandle.
    # Otherwise, you will have too many open operating system filehandles
    # if you run this function many times (more than ulimit -n allows, that is)
    empirical_list = myEmpiricalList['arr_0']
    myEmpiricalList.close()

    return empirical_list


def pvalue_from_empirical_list(tandemrepeat, score_type='phylo', score_value=None, empirical=[]):
    """Calculates the p-Value of a score_type for the given `tandemrepeat`.


    The p-Value is the number of scores for comparable tandem repeats, i.e. of
    same repeat unit length and repeat unit copy number that are as good or
    better.

    Args:
        tandemrepeat (Repeat): An instance of the Repeat class.
        score_type: The model of repeat evolution used. E.g. "phylo".
        score_value (float): The value of the score of the Repeat instance.
        empirical: The null distribution of the score for repeats of the same
                    type.

    Returns:
        pvalue: A float.

    """
    if score_value is None:
        score_value = tandemrepeat.score(score_type)

    if len(empirical) == 0:
        empirical = empirical_list(tandemrepeat.l_effective, tandemrepeat.n,
                                  tandemrepeat.sequence_type, score_type)

    # A smaller score is a better score for all scores in this list.
    if score_type in ['entropy']:
        if score_value == -1:
            return 1
        else:
            return bisect.bisect_right(empirical, score_value) / len(empirical)
    else:
        return 1 - bisect.bisect_left(empirical, score_value) / len(empirical)


# ####################### pSim & parsimony: read in PDF #######################

def column_pdf(n, score_type='psim', sequence_type='AA'):
    """Load and return the probability density function of the score on random
        ``sequence_type`` data of length ``n``.

    Load and return the probability density function of the score on random
    ``sequence_type`` data of length ``n``.
    This method handels scores with distributions independent of ``l``. This includes
    all parsimony scores, and excludes model based scores.

    The order of the values in the probability density function is WORST FIRST: The first
    values describe the probability for "bad" scores, the last values the probabilities
    for "top" scores.

    Currently, for all models values are available until

        * ``n_max`` = 150

    For ``n`` <=  ``n_max`` all values are available. If no distribution is available for
    a given ``n`` the closest values for ``n`` is used.

    For values of ``n`` above ``n_max``, the ``n_max`` pdf is returned.

    Args:
        n (int): The number of repeat units.
        score_type: The heuristic model used. E.g. "psim".
        sequence_type (str): The type of the sequence: either "AA" or "DNA".

    Returns:
        A numpy list of length 10.000

    .. todo:: Check return type.
    .. todo:: Check how this method concurs with the ``empirical_list`` method.
    """

    # Currently, the pdfs are only save up to n_max = 150
    n_max = 150
    with open(join(path_score, sequence_type, score_type, str(min(n, n_max)) + '.txt'), 'r') as pdf_file:
        pdf = [i for i in csv.reader(pdf_file, dialect='excel-tab')]
        pdf = [np.double(i[0]) for i in pdf[1:]]
        return np.array(pdf)

# ############################# pSim & parsimony ##############################


def calculate_repeat_structure(tandemrepeat):
    """ Calculate the number of columns with a certain number of gaps for each column in
    a ``Repeat`` instance.

    You can use a different null distribution for each column of different length for both
    the parsimony and the pSim score. Therefore, calculate the structure of the TR before
    you calculate the null distribution of tandem repeat scores of tandem repeats with
    the same gap distribution.

    Args:
        n (int): The number of repeat units.
        score: The heuristic model used. E.g. "psim".
        sequence_type (str): The type of the sequence: either "AA" or "DNA".

    Returns:
        l_repeat_structure (list)
        n_repeat_structure (int)

    .. todo:: Include exact definition of returned values.
    """

    repeat_structure = [
        len(column.replace("-", ""))
        for column in tandemrepeat.msaTD if (2 * column.count('-') < len(column))
    ]
    n_repeat_structure = list(set(repeat_structure))
    l_repeat_structure = [
        repeat_structure.count(i) for i in n_repeat_structure
    ]
    return l_repeat_structure, n_repeat_structure


# #################### DERIVE p-Value distributions ###########################

def calc_pvalues(repeats, result_file_path, file_name,
        scoretypes=['phylo', 'phylo_gap'], gappy_data=False):
    """ Create and save a null distribution for ``repeats`` scores for p-Value calculation.

    You can use a different null distribution for each column of different length for both
    the parsimony and the pSim score. Therefore, calculate the structure of the TR before
    you calculate the null distribution of tandem repeat scores of tandem repeats with
    the same gap distribution.

    Args:
        repeats (list of Repeat): A list of ``Repeat`` instances.
        result_file_path (str): Path to result folder.
        file_name (str): Name of result file.
        scoretypes (list of str):  List of scores. E.g. ['phylo','phylo_gap'].
        gappy_data (bool): True if any of ``repeats`` contain gaps, else False.

    """

    # If the repeats are not gappy, than you can use always the same distribution of scores
    # on random data to calculate the p-Value. Otherwise, there might be deletion columns
    # and the distribution should be loaded each time
    # 'pvalue_from_empirical_list' is called.
    empirical_list = []

    for i_score_type in scoretypes:
        if 'parsimony' == i_score_type:
            test_statistic = [
                pvalue_pars(i_repeat) if i_repeat.l_effective != 0 else 1 for i_repeat in repeats]
        elif 'pSim' == i_score_type:
            test_statistic = [
                pvalue_psim(i_repeat) if i_repeat.l_effective != 0 else 1 for i_repeat in repeats]
        else:
            if not gappy_data:
                empirical_list = empirical_list(
                    repeats[0].l_effective,
                    repeats[0].n,
                    repeats[0].sequence_type,
                    i_score_type)
            test_statistic = [
                pvalue_from_empirical_list(
                    i_repeat,
                    i_score_type,
                    empirical=empirical_list) if i_repeat.l_effective != 0 else 1 for i_repeat in repeats]
        score_path = os.path.join(result_file_path, i_score_type)
        if not os.path.isdir(score_path):
            os.makedirs(score_path)
        print(os.path.join(score_path, file_name))
        np.savez(os.path.join(score_path, file_name), np.sort(test_statistic))


def d_average_multinom(l, n, sequence_type, score_type):
    """ Helper function for the analytic calculation of parsimony or pSim scores.

    Called by ``dAverageMultipleMaxMultinom`` or ``d_average_multiple_pars_multinom``.

    Return the ``l``-times self-convoluted pdf of the score on random sequence_type data
    of length ``n``. The order of the pdf is kept from ``column_pdf()``.

    For details see:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        l (int): The length of the repeat unit.
        n (int): The number of repeat units.
        sequence_type (str): The type of the sequence: either "AA" or "DNA".
        score_type: The model of repeat evolution used. Either 'psim' or 'parsimony'.

    Returns:
        (description missing)

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html

    .. todo:: Describe return value.
    """

    complete_pdf = column_pdf(n=n, score_type=score_type, sequence_type=sequence_type)
    if l != 1:
        single_column_pdf = complete_pdf
        for i in range(2, l + 1):
            complete_pdf = sp.convolve(complete_pdf, single_column_pdf)
    return complete_pdf

# ####################### pSim ################################################


# python 2: precision must be float
def dAverageMultipleMaxMultinom(tandemrepeat, precision=10000.):
    """ Calculate null distribution for e.g. the pSim score for repeats of type ``tandemrepeat``
        as a probability density function.

    Analytically calculate the p-Value distribution for the pSim score for repeats of
    type ``tandemrepeat``. The derivation is described in:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        tandemrepeat (Repeat): A ``Repeat`` instance.
        precision (float): The precision of the returned probability density function in
        terms of the length of the resulting list.

    Returns:
        p (list of float): cumulated probabilities from 0 to 1.
        unnamed (list of float): scores corresponding to the probabilities in ``p``.

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html
    """

    l_repeat_structure, n_repeat_structure = calculate_repeat_structure(tandemrepeat)

    p = d_average_multinom(
        l_repeat_structure[0], n_repeat_structure[0],
        tandemrepeat.sequence_type, score='psim'
    )
    val = np.array(
        (np.r_[l_repeat_structure[0]:(l_repeat_structure[0] * n_repeat_structure[0] + 1.)]
         * (precision / n_repeat_structure[0])).round(),
        dtype='uint32'
    )
    if len(l_repeat_structure) == 1:
        # return best values first. for pSim, best = 1, worst = 0
        return p[::-1], val[::-1] / (precision * sum(l_repeat_structure))
    else:
        for i in range(1, len(l_repeat_structure)):
            p = np.outer(p,
                         d_average_multinom(
                             l_repeat_structure[i], n_repeat_structure[i],
                             tandemrepeat.sequence_type, score='psim')
                         ).ravel()

            val = (val[:, np.newaxis] +
                   np.array(
                (np.r_[l_repeat_structure[i]:(l_repeat_structure[i] * n_repeat_structure[i] + 1.)]
                 * (precision / n_repeat_structure[i])
                 ).round(),
                dtype='uint32')
            ).ravel()

            x = np.bincount(val, weights=p)
            val = np.unique(val)
            p = x[val]
        # return best values first. for pSim, best = 1, worst = 0
        return p[::-1], (val / (precision * sum(l_repeat_structure)))[::-1]


def pvalue_psim(tandemrepeat):
    """ Calculate the p-Value of the pSim score for a Repeat.

    Retrieve the probability density function for repeats of the same type as ``tandemrepeat``.
    Then, calculate the p-Value given this probability density function, and the
    pSim score of ``tandemrepeat``.

    Args:
        tandemrepeat (Repeat): A ``Repeat`` instance.

    Returns:
        p-Value (float)

    .. todo:: Check the method's behaviour if ``tandemrepeat`` s parsimony score has not been
        calculated before.
    .. todo:: Check exception: if pdf == False: return 1.
    .. todo:: Describe ``precision``.
    """

    precision = 10000.

    pdf = dAverageMultipleMaxMultinom(tandemrepeat, precision)
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(tandemrepeat.score('pSim') == pdf[1])[0]
    index = np.where(np.abs(tandemrepeat.score('pSim') - pdf[1]) <= 1. / precision)[0]
    if len(index) == 1:  # standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(tandemrepeat.score('pSim') > pdf[1])[0]
    if len(indices) >= 1:
        # if pars is not exactly included in list, give back mean of value above
        # and value below (should only occur due to numerical imprecision)
        a = min(indices)
        return((cumsumPDF[a] + cumsumPDF[a - 1]) / 2)
    # This should only occur if score('pSim') is really bad or really good
    else:
        return 1

# ################## parsimony ################################################


# python 2: precision must be float
def d_average_multiple_pars_multinom(tandemrepeat, precision=10000.):
    """ Calculate null distribution for the parsimony score for repeats of type ``tandemrepeat``
        as a probability density function.

    Analytically calculate the p-Value distribution for the parsimony score for repeats of
    type ``tandemrepeat``. The derivation is described in:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        tandemrepeat (Repeat): A ``Repeat`` instance.
        precision (float): The precision of the returned probability density function in
        terms of the length of the resulting list.

    Returns:
        p (list of float): cumulated probabilities from 0 to 1.
        unnamed (list of float): scores corresponding to the probabilities in ``p``.

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html
    """

    l_repeat_structure, n_repeat_structure = calculate_repeat_structure(tandemrepeat)

    p = d_average_multinom(
        l_repeat_structure[0], n_repeat_structure[0],
        tandemrepeat.sequence_type, score='parsimony'
    )
    val = np.array(
        (np.r_[(l_repeat_structure[0] * (n_repeat_structure[0] - 1)):-1:-1.]
         * (precision / (n_repeat_structure[0] - 1.))
         ).round(),
        dtype='uint32'
    )
    if len(l_repeat_structure) == 1:
        # return best values first. for parsimony, best = 0, worst = 1
        return p[::-1], (val / (precision * sum(l_repeat_structure)))[::-1]
    else:
        for i in range(1, len(l_repeat_structure)):
            try:
                p = np.outer(p,
                             d_average_multinom(
                                 l_repeat_structure[i],
                                 n_repeat_structure[i],
                                 tandemrepeat.sequence_type, score='parsimony')
                             ).ravel()

                val = (val[:, np.newaxis] + np.array(
                    (np.r_[(l_repeat_structure[i] * (n_repeat_structure[i] - 1)):-1:-1.]
                     * (precision / (n_repeat_structure[i] - 1.))
                     ).round(),
                    dtype='uint32')
                ).ravel()

                x = np.bincount(val, weights=p)
                val = np.unique(val)
                p = x[val]
            except:
                LOG.warning(
                    "Failed on: " + str(l_repeat_structure) + " " +
                    str(n_repeat_structure)
                )
                return False
        # return best values first. for parsimony, best = 0, worst = 1.
        # Here, the  np.bincount() and np.unique() funs did the reordering.
        return p, val / (precision * sum(l_repeat_structure))


def pvalue_pars(tandemrepeat):
    """ Calculate the p-Value of the parsimony score for a Repeat.

    Retrieve the probability density function for repeats of the same type as ``tandemrepeat``.
    Then, calculate the p-Value given this probability density function, and the
    parsimony score of ``tandemrepeat``.

    Args:
        tandemrepeat (Repeat): A ``Repeat`` instance.

    Returns:
        p-Value (float)

    .. todo:: Check the method's behaviour if ``tandemrepeat`` s parsimony score has not been
        calculated before.
    .. todo:: Check exception: if pdf == False: return 1.
    """

    precision = 10000.

    pdf = d_average_multiple_pars_multinom(tandemrepeat, precision)
    # Check the following three lines:
    if not pdf:
        return 1
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(pdf[1] == tandemrepeat.score('parsimony'))[0]
    index = np.where(
        np.abs(
            tandemrepeat.score('parsimony') -
            pdf[1]) <= 1. /
        precision)[0]
    if len(index) == 1:  # standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(tandemrepeat.score('parsimony') < pdf[1])[0]
    if len(indices) >= 1:
        # if pars is not exactly included in list, give back mean of value
        # above and value below (should only occur due to numerical
        # imprecision)
        a = min(indices)
        return (cumsumPDF[a] + cumsumPDF[a - 1]) / 2
    else:  # This should only occur if score('parsimony') is really bad
        return 1


# #################### gap penalty ############################################

def gap_penalty(tandemrepeat, mu):
    """ Calculate the gap penalty for a ``Repeat`` given mutation rate ``mu``.

    Args:
        tandemrepeat (Repeat): A ``Repeat`` instance.
        mu (float): The mutation rate.

    Returns:
        Gap penalty (float)

    .. todo:: Define ``mu`` more precisely.
    .. todo:: Is this function called from anywhere? In case, consider refactoring.
    """

    return ((mu / (1 - mu)) ** tandemrepeat.n_gap_structure) * tandemrepeat.p_gap_structure
