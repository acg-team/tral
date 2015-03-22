# (C) 2014 Elke Schaper

import os
import bisect
from os.path import join
import csv
import logging
import numpy as np
import scipy as sp
import scipy.stats
import scipy.special

from tral.paths import DATA_DIR

log = logging.getLogger(__name__)

path_score = join(DATA_DIR, 'pvalue')

# ######################### REPEAT SCORE P-VALUE CALCULATION FUNCTIONS ########

# #################### phylo & entropy ########################################


def empiricalList(l, n, sequence_type='AA', score='phylo'):
    """Load and return a numpy list with 10,000 empirical values of the user defined
        distribution from an external file.


    If no distribution is available for a given ``l`` and ``n``, the closest values for
    ``l`` and ``n`` are chosen instead.

    For the standard model score "phylo", currently values are available until

        * ``lMax`` = 99
        * ``nMax`` = 49
        * ``total_repeat_length_max`` = 1000

    For ``l`` <=  ``lMax`` and ``n`` <=  ``nMax`` and ``n*l`` <=  ``total_repeat_length_max``
    all values are available.

    For other values, we use the closest distribution, assuming that values change little
    (heuristic assumption!).

    Args:
        l (int): The length of the repeat unit.
        n (int): The number of repeat units.
        sequence_type (str): The type of the sequence: either "AA" or "DNA".
        score: The model of repeat evolution used. E.g. "phylo".

    Returns:
        A numpy list of length 10.000


    .. todo:: Define "phylo" model.
    .. todo:: Perhaps add distributions beyond lMax, nMax.
    """

    lMax = 99
    nMax = 49
    if score == 'entropy':
        lMax = 20
        nMax = 5

    total_repeat_length_max = 1000
    l = min(l, lMax)
    n = min(n, nMax)
    if l * n > total_repeat_length_max:
        logging.debug(
            "l: %d and n: %d are bigger than the total_repeat_length_max: %d" %
            (l, n, total_repeat_length_max))
        # Dirty little hack: as we do not have data for this pair of l and n,
        # apply nearest data file.
        while l * n > total_repeat_length_max:
            if l > n:
                l -= 1
            else:
                n -= 1

    file = join(
        path_score,
        sequence_type,
        score,
        str(l) +
        '_' +
        str(n) +
        '.npz')
    if not os.path.isfile(file):
        raise ValueError("complete pdf file %s does not exist!" % file)

    myEmpiricalList = np.load(file)
    # It is absolutely necessary to close the numpy filehandle.
    # Otherwise, you will have too many open operating system filehandles
    # if you run this function many times (more than ulimit -n allows, that is)
    empiricalList = myEmpiricalList['arr_0']
    myEmpiricalList.close()

    return empiricalList


def pvalueFromEmpiricialList(myTR, score='phylo', myScore=None, empirical=[]):
    """Calculates the p-Value of a score for the given myTR.


    The p-Value is the number of scores for comparable tandem repeats, i.e. of same
    repeat unit length and repeat unit copy number that are as good or better.

    Args:
        myTR (Repeat): An instance of the Repeat class.
        score: The model of repeat evolution used. E.g. "phylo".
        myScore (float): The value of the score of the Repeat instance.
        empirical: The null distribution of the score for repeats of the same type.

    Returns:
        pvalue: A float.

    """
    if myScore is None:
        myScore = myTR.score(score)

    if len(empirical) == 0:
        empirical = empiricalList(myTR.l_effective, myTR.n, myTR.sequence_type, score)

    # A smaller score is a better score for all scores in this list.
    if score in ['entropy']:
        if myScore == -1:
            return 1
        else:
            return bisect.bisect_right(empirical, myScore) / len(empirical)
    else:
        return 1 - bisect.bisect_left(empirical, myScore) / len(empirical)


########################### pSim & parsimony: read in PDF ################

def columnPDF(n, score='psim', sequence_type='AA'):
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

        * ``nMax`` = 150

    For ``n`` <=  ``nMax`` all values are available. If no distribution is available for
    a given ``n`` the closest values for ``n`` is used.

    For values of ``n`` above ``nMax``, the ``nMax`` pdf is returned.

    Args:
        n (int): The number of repeat units.
        score: The heuristic model used. E.g. "psim".
        sequence_type (str): The type of the sequence: either "AA" or "DNA".

    Returns:
        A numpy list of length 10.000

    .. todo:: Check return type.
    .. todo:: Check how this method concurs with the ``empiricalList`` method.
    """

    # Currently, the pdfs are only save up to nMax = 150
    nMax = 150
    with open(join(path_score, sequence_type, score, str(min(n, nMax)) + '.txt'), 'r') as pdf_file:
        pdf = [i for i in csv.reader(pdf_file, dialect='excel-tab')]
        pdf = [np.double(i[0]) for i in pdf[1:]]
        return np.array(pdf)

################################### pSim & parsimony #####################


def calculate_repeat_structure(myTR):
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
        lRepeatStructure (list)
        nRepeatStructure (int)

    .. todo:: Include exact definition of returned values.
    """

    repeatStructure = [
        len(column.replace("-", ""))
        for column in myTR.msaTD if (2 * column.count('-') < len(column))
    ]
    nRepeatStructure = list(set(repeatStructure))
    lRepeatStructure = [
        repeatStructure.count(i) for i in nRepeatStructure
    ]
    return lRepeatStructure, nRepeatStructure


############################ DERIVE p-Value distributions ################

def calc_pvalues(
        repeats,
        resultFilePath,
        fileName,
        scoreslist=[
            'phylo',
            'phylo_gap'],
        gappy_data=False):
    """ Create and save a null distribution for ``repeats`` scores for p-Value calculation.

    You can use a different null distribution for each column of different length for both
    the parsimony and the pSim score. Therefore, calculate the structure of the TR before
    you calculate the null distribution of tandem repeat scores of tandem repeats with
    the same gap distribution.

    Args:
        repeats (list of Repeat): A list of ``Repeat`` instances.
        resultFilePath (str): Path to result folder.
        fileName (str): Name of result file.
        scoreslist (list of str):  List of scores. E.g. ['phylo','phylo_gap'].
        gappy_data (bool): True if any of ``repeats`` contain gaps, else False.

    """

    # If the repeats are not gappy, than you can use always the same distribution of scores
    # on random data to calculate the p-Value. Otherwise, there might be deletion columns
    # and the distribution should be loaded each time
    # 'pvalueFromEmpiricialList' is called.
    empirical_list = []

    for iScore in scoreslist:
        if 'parsimony' == iScore:
            testStatistic = [
                pvaluePars(iRepeat) if iRepeat.l_effective != 0 else 1 for iRepeat in repeats]
        elif 'pSim' == iScore:
            testStatistic = [
                pvaluePSim(iRepeat) if iRepeat.l_effective != 0 else 1 for iRepeat in repeats]
        else:
            if not gappy_data:
                empirical_list = empiricalList(
                    repeats[0].l_effective,
                    repeats[0].n,
                    repeats[0].sequence_type,
                    iScore)
            testStatistic = [
                pvalueFromEmpiricialList(
                    iRepeat,
                    iScore,
                    empirical=empirical_list) if iRepeat.l_effective != 0 else 1 for iRepeat in repeats]
        score_path = os.path.join(resultFilePath, iScore)
        if not os.path.isdir(score_path):
            os.makedirs(score_path)
        print(os.path.join(score_path, fileName))
        np.savez(os.path.join(score_path, fileName), np.sort(testStatistic))


def dAverageMultinom(l, n, sequence_type, score):
    """ Helper function for the analytic calculation of parsimony or pSim scores.

    Called by ``dAverageMultipleMaxMultinom`` or ``dAverageMultipleParsMultinom``.

    Return the ``l``-times self-convoluted pdf of the score on random sequence_type data
    of length ``n``. The order of the pdf is kept from ``columnPDF()``.

    For details see:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        l (int): The length of the repeat unit.
        n (int): The number of repeat units.
        sequence_type (str): The type of the sequence: either "AA" or "DNA".
        score: The model of repeat evolution used. Either 'psim' or 'parsimony'.

    Returns:
        (description missing)

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html

    .. todo:: Describe return value.
    """

    completePDF = columnPDF(n=n, score=score, sequence_type=sequence_type)
    if l != 1:
        singleColumnPDF = completePDF
        for i in range(2, l + 1):
            completePDF = sp.convolve(completePDF, singleColumnPDF)
    return completePDF

####################################### pSim #############################


# python 2: precision must be float
def dAverageMultipleMaxMultinom(myTR, precision=10000.):
    """ Calculate null distribution for e.g. the pSim score for repeats of type ``myTR``
        as a probability density function.

    Analytically calculate the p-Value distribution for the pSim score for repeats of
    type ``myTR``. The derivation is described in:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        myTR (Repeat): A ``Repeat`` instance.
        precision (float): The precision of the returned probability density function in
        terms of the length of the resulting list.

    Returns:
        p (list of float): cumulated probabilities from 0 to 1.
        unnamed (list of float): scores corresponding to the probabilities in ``p``.

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html
    """

    lRepeatStructure, nRepeatStructure = calculate_repeat_structure(myTR)

    p = dAverageMultinom(
        lRepeatStructure[0], nRepeatStructure[0],
        myTR.sequence_type, score='psim'
    )
    val = np.array(
        (np.r_[lRepeatStructure[0]:(lRepeatStructure[0] * nRepeatStructure[0] + 1.)]
         * (precision / nRepeatStructure[0])).round(),
        dtype='uint32'
    )
    if len(lRepeatStructure) == 1:
        # return best values first. for pSim, best = 1, worst = 0
        return p[::-1], val[::-1] / (precision * sum(lRepeatStructure))
    else:
        for i in range(1, len(lRepeatStructure)):
            p = np.outer(p,
                         dAverageMultinom(
                             lRepeatStructure[i], nRepeatStructure[i],
                             myTR.sequence_type, score='psim')
                         ).ravel()

            val = (val[:, np.newaxis] +
                   np.array(
                (np.r_[lRepeatStructure[i]:(lRepeatStructure[i] * nRepeatStructure[i] + 1.)]
                 * (precision / nRepeatStructure[i])
                 ).round(),
                dtype='uint32')
            ).ravel()

            x = np.bincount(val, weights=p)
            val = np.unique(val)
            p = x[val]
        # return best values first. for pSim, best = 1, worst = 0
        return p[::-1], (val / (precision * sum(lRepeatStructure)))[::-1]


def pvaluePSim(myTR):
    """ Calculate the p-Value of the pSim score for a Repeat.

    Retrieve the probability density function for repeats of the same type as ``myTR``.
    Then, calculate the p-Value given this probability density function, and the
    pSim score of ``myTR``.

    Args:
        myTR (Repeat): A ``Repeat`` instance.

    Returns:
        p-Value (float)

    .. todo:: Check the method's behaviour if ``myTR`` s parsimony score has not been
        calculated before.
    .. todo:: Check exception: if pdf == False: return 1.
    .. todo:: Describe ``precision``.
    """

    precision = 10000.

    pdf = dAverageMultipleMaxMultinom(myTR, precision)
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(myTR.score('pSim') == pdf[1])[0]
    index = np.where(np.abs(myTR.score('pSim') - pdf[1]) <= 1. / precision)[0]
    if len(index) == 1:  # standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(myTR.score('pSim') > pdf[1])[0]
    if len(indices) >= 1:
        # if pars is not exactly included in list, give back mean of value above
        # and value below (should only occur due to numerical imprecision)
        a = min(indices)
        return((cumsumPDF[a] + cumsumPDF[a - 1]) / 2)
    # This should only occur if score('pSim') is really bad or really good
    else:
        return 1

#################################### parsimony ###########################


# python 2: precision must be float
def dAverageMultipleParsMultinom(myTR, precision=10000.):
    """ Calculate null distribution for the parsimony score for repeats of type ``myTR``
        as a probability density function.

    Analytically calculate the p-Value distribution for the parsimony score for repeats of
    type ``myTR``. The derivation is described in:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
        myTR (Repeat): A ``Repeat`` instance.
        precision (float): The precision of the returned probability density function in
        terms of the length of the resulting list.

    Returns:
        p (list of float): cumulated probabilities from 0 to 1.
        unnamed (list of float): scores corresponding to the probabilities in ``p``.

    .. warning:: if precision higher than max(uint32) use uint64 instead.
            CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html
    """

    lRepeatStructure, nRepeatStructure = calculate_repeat_structure(myTR)

    p = dAverageMultinom(
        lRepeatStructure[0], nRepeatStructure[0],
        myTR.sequence_type, score='parsimony'
    )
    val = np.array(
        (np.r_[(lRepeatStructure[0] * (nRepeatStructure[0] - 1)):-1:-1.]
         * (precision / (nRepeatStructure[0] - 1.))
         ).round(),
        dtype='uint32'
    )
    if len(lRepeatStructure) == 1:
        # return best values first. for parsimony, best = 0, worst = 1
        return p[::-1], (val / (precision * sum(lRepeatStructure)))[::-1]
    else:
        for i in range(1, len(lRepeatStructure)):
            try:
                p = np.outer(p,
                             dAverageMultinom(
                                 lRepeatStructure[i],
                                 nRepeatStructure[i],
                                 myTR.sequence_type, score='parsimony')
                             ).ravel()

                val = (val[:, np.newaxis] + np.array(
                    (np.r_[(lRepeatStructure[i] * (nRepeatStructure[i] - 1)):-1:-1.]
                     * (precision / (nRepeatStructure[i] - 1.))
                     ).round(),
                    dtype='uint32')
                ).ravel()

                x = np.bincount(val, weights=p)
                val = np.unique(val)
                p = x[val]
            except:
                logging.warning(
                    "Failed on: " + str(lRepeatStructure) + " " +
                    str(nRepeatStructure)
                )
                return False
        # return best values first. for parsimony, best = 0, worst = 1.
        # Here, the  np.bincount() and np.unique() funs did the reordering.
        return p, val / (precision * sum(lRepeatStructure))


def pvaluePars(myTR):
    """ Calculate the p-Value of the parsimony score for a Repeat.

    Retrieve the probability density function for repeats of the same type as ``myTR``.
    Then, calculate the p-Value given this probability density function, and the
    parsimony score of ``myTR``.

    Args:
        myTR (Repeat): A ``Repeat`` instance.

    Returns:
        p-Value (float)

    .. todo:: Check the method's behaviour if ``myTR`` s parsimony score has not been
        calculated before.
    .. todo:: Check exception: if pdf == False: return 1.
    """

    precision = 10000.

    pdf = dAverageMultipleParsMultinom(myTR, precision)
    # Check the following three lines:
    if not pdf:
        return 1
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(pdf[1] == myTR.score('parsimony'))[0]
    index = np.where(
        np.abs(
            myTR.score('parsimony') -
            pdf[1]) <= 1. /
        precision)[0]
    if len(index) == 1:  # standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(myTR.score('parsimony') < pdf[1])[0]
    if len(indices) >= 1:
        # if pars is not exactly included in list, give back mean of value
        # above and value below (should only occur due to numerical
        # imprecision)
        a = min(indices)
        return (cumsumPDF[a] + cumsumPDF[a - 1]) / 2
    else:  # This should only occur if score('parsimony') is really bad
        return 1


####################################### gap penalty ######################

def gapPenalty(myTR, mu):
    """ Calculate the gap penalty for a ``Repeat`` given mutation rate ``mu``.

    Args:
        myTR (Repeat): A ``Repeat`` instance.
        mu (float): The mutation rate.

    Returns:
        Gap penalty (float)

    .. todo:: Define ``mu`` more precisely.
    .. todo:: Is this function called from anywhere? In case, consider refactoring.
    """

    return ((mu / (1 - mu)) ** myTR.nGapStructure) * myTR.pGapStructure
