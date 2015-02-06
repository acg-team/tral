# (C) 2011 Alexander Korsunsky
# (C) 2011-2014 Elke Schaper

import collections
import logging
import math
import numpy as np
import os
from os.path import join
import re
import scipy as sp
import scipy.stats, scipy.special, scipy.linalg

from tral import configuration
from tral.paths import *

log = logging.getLogger(__name__)

c = configuration.Configuration.Instance()
config_general = c.config
config = config_general["repeat_score"]

########################## REPEAT SCORE CALCULATION FUNCTIONS ############################

############################## 1. SIMILARITY SCORES ######################################

def load_equilibrium_freq(filename):
    """ Load equilibrium frequencies from a substitution matrix file.

    Load equilibrium frequencies from a substitution matrix file. Note that with minor
     changes, the whole amino acid substitution matrix could be returned.

    Args:
        filename (str): Path to substitution matrix file, E.g. LG.dat.

    Returns:
        (Dict of str,float): Dictionary of one letter amino acid codes and float value
    """
    patstr_freqline = r"[+-]?\d+(?:\.\d+)?\s*"
    with open(filename,"r") as f:
        # Create empty 20x20 matrix
        substitution_matrix = [[0.0] * 20] * 20
        row_count = 0
        for row in f:
            #print("row %d: " % row_count, row)
            row = row.strip()
            pat = r"\s*" + patstr_freqline * (row_count+1)
            match = re.findall(pat, row)
            if not match:
                continue
            # store read lines into matrix
            substitution_matrix[row_count][0:row_count+1] = map(float, row.split())
            row_count = row_count+1
            # Read 19 lines of substitution matrix
            if row_count == 19:
                break
        else:
            # EOF before reading the frequency line
            raise ValueError("Input file is malformed! Did not find enough lines for the substitution matrix.")
        frequencies = []
        for row in f:
            matches = re.findall(patstr_freqline, row)
            frequencies.extend(map(float, matches))
            if len(frequencies) >= 20:
                break
        else:
            # EOF before finishing the frequency line
            raise ValueError("Input file is malformed! Did not find enough numbers for equilibrium frequencies.")
        substitution_matrix[19][0:20] = frequencies[0:20]
    #return substitution_matrix
    return { aa : freq for aa, freq in
        zip(
            ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"],
            substitution_matrix[19]
        )
    }


def meanSimilarity(repeat, measureOfSimilarity, ignoreGaps = config['indel'].as_bool('ignore_gaps')):
    """ Calculate the mean similarity of all columns in the multiple sequence alignment
    in ``repeat`` according to a ``measureOfSimilarity``.

    Args:
        repeat (Repeat): An instance of the Repeat class.
        measureOfSimilarity (str): Either "entropy", "parsimony", or "pSim".
        ignoreGaps (boolean): Are gaps, "-" treated as chars (False) or ignored (True).

    Returns:
        float: A similarity measure for the repeat units in ``repeat``.

    ..  todo:: Is ``rep_and_meas`` used at current.
    """

    if not ignoreGaps:
        return sum(map(measureOfSimilarity, repeat.msaTD)) / float(repeat.l)
    else:
        def rep_and_meas(column):
            return measureOfSimilarity(column.replace("-", ""))

        return sum(map(
                lambda column: measureOfSimilarity(column.replace("-", "")),
                repeat.msaTD)) / repeat.lD  ## python 2: integer division

def entropy(column):
    """ Calculate the entropy of a list.

     Calculate the entropy of list. Assume that the frequency of elements in the list
     reflects a probability density estimator.

    For an exact description see:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
      column (list of str): List of elements that can be ordered in a set, e.g. str or int.

    Returns:
        float: The entropy of ``column``.
    """

    estimatedP = [(column.count(a)/float(len(column))) for a in set(column)]
    return sum([-p*np.log2(p) for p in estimatedP if p > 0])

def parsimony(column):
    """ Calculate the parsimony score of a list.

    Calculate the parsimony score on ``column``:  Count the number of changes in a
    ``column`` of a multiple sequence alignment and divide by the maximally possible
    number of changes.

    For an exact description see:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
      column (list of str): A ``list`` of two or more elements that can be ordered in a ``set``, e.g.
       ``str`` or ``int``

    Returns:
        float: The parsimony score of ``column``.

    ..  todo:: At which step should we assert that the length of ``column`` > 1?
    """

    if len(column) > 1:
        # Here, 0 means no changes and 1 only changes # python 2: integer division
        return (len(set(column))-1)/(len(column)-1)
    else:
        logging.warning("Repeats with n = 1 are not repeats")
        #raise AssertionError("Repeats with n = 1 are not repeats")
        return 1

def pSim(column):
    """ Calculate the pSim score of a list.

    Calculate the pSim score on ``column``:  Count the frequency of the most frequent
    element in a ``column`` of a multiple sequence alignment.

    For an exact description see:

    Schaper, E., Kajava, A., Hauser, A., & Anisimova, M. Repeat or not repeat?
    --Statistical validation of tandem repeat prediction in genomic sequences.
    Nucleic Acids Research (2012).

    Args:
      column (list of str): A list of elements.

    Returns:
        float: The pSim score of ``column``.
    """

    return collections.Counter(column).most_common(1)[0][1]/float(len(column))



#################### 2. MAXIMUM LIKELIHOOD CALCULATIONS ##################################


def optimisation(function, args, start_min = config['optimisation'].as_float('start_min'),
                start_max = config['optimisation'].as_float('start_max'),
                nIteration = config['optimisation'].as_int('nIteration')):

    """Binary one-dimensional optimisation of ``function``.

    Perform a binary one dimensional optimisation of parameter ``t`` on ``function`` with
    the additional list of arguments/parameters ``args``. Start on a given range of start,
    and iterate ``nIteration`` times. Return the maximum of the function, together with
    the maximum parameter ``t`` as a tuple.

    In addition to ``start_min`` and ``start_max``, additional hard boundaries
    on the optimisation parameter ``t`` are hard-coded: t E [0.0000000001, 100]

    Args:
      function (method): Method with input (float, ``args``) and output float.
      args (list of items): List of arguments for ``function``.

    Kwargs:
        start_min (float): Minimum value of function parameter ``t``.
        start_max (float): Maximum value of function parameter ``t``.
        nIteration (int): Number of iterations until optimisation results are returned.

    Returns:
        tuple (Maximum of function: float, Corresponding function parameter: float)

      ..  todo:: Check why ``regionBounded`` is hardcoded at current.
      ..  todo:: Check why ``stepsize`` is hardcoded at current.
      ..  todo:: Generalise hard-coded optimisation parameter boundaries.
    """

    t = [[start_min], [(start_min+start_max)/2], [start_max]]
    stepsize = 2
    for iT in t:
        iT.append(function(iT[0], *args))
    regionBounded = False
    count = 0
    while count < nIteration:
        count += 1
        if not regionBounded:
            # First entry is best - shift trial out zone to left
            if t[0][1] > t[1][1]:
                t[2] = t[1][:]
                t[1] = t[0][:]
                t[0] = [  max(0.0000000001,t[0][0] - stepsize)   ]
                t[0].append(function(t[0][0], *args))
            # Last entry is best - shift trial out zone to right
            elif t[2][1] > t[1][1]:
                t[0] = t[1][:]
                t[1] = t[2][:]
                t[2] = [min(100,t[2][0] + stepsize)]
                t[2].append(function(t[2][0], *args))
            else:
                regionBounded = True
        if regionBounded:
            # calc halfPoints
            halfPoints = [[(t[0][0] + t[1][0])/2], [(t[1][0] + t[2][0])/2]]
            for iT in halfPoints:
                iT.append(function(iT[0], *args))
            # First halfPoint is best - shift it to be the next midPoint
            if halfPoints[0][1] > t[1][1]:
                t[2] = t[1][:]
                t[1] = halfPoints[0][:]
            # Last halfPoint is best - shift it to be the next midPoint
            elif halfPoints[1][1] > t[1][1]:
                t[0] = t[1][:]
                t[1] = halfPoints[1][:]
            # Former midPoint is still the best - adjust the boundaries
            else:
                t[0] = halfPoints[0][:]
                t[2] = halfPoints[1][:]
    return (t[1][0], t[1][1])

def loadModel(evolution_model = 'lg'):
    ## Prepare the model of repeat evolution
    if evolution_model == 'lg':
        # 1. Substitutions are modeled by the LG matrix
        aa,eqFreq,subsMatrix = LG()

    else:
        aa,eqFreq,subsMatrix = K80()

    equilibrium_freq = {i:j for i,j in zip(aa,eqFreq)}
    alphabet = { b:a for a,b in enumerate(aa)}
    length = len(aa)
    eqFreq = np.array(eqFreq)
    eqFreq = eqFreq/sum(eqFreq)
    S = np.array(subsMatrix)
    PI = np.diag(eqFreq)
    Q = np.dot(S,PI)
    for i in range(length):
        Q[i,i] = - sum(Q[i][:i]) - sum(Q[i][i+1:])
    diagonal = [Q[i,i] for i in range(length)]
    Q *= - sum( (i*j for i,j in zip(diagonal, eqFreq) ))
    return (Q,eqFreq,alphabet)

def TN93(alpha_1 = config['TN93']['alpha_1'], alpha_2 = config['TN93']['alpha_2'], beta = config['TN93']['beta']):
    '''TN93'''
    #4*4 matrix with α1 = 0.3, α2 = 0.4, β = 0.7 according to TN93
    Q = [[-(alpha_1 + 2*beta), alpha_1, beta, beta],[alpha_1, -(alpha_1 + 2*beta), beta, beta], \
    [beta, beta, -(alpha_2 + 2*beta), alpha_2], [beta, beta, alpha_2, -(alpha_2 + 2*beta)]]

    return(['T', 'C', 'A','G'],[0.25, 0.25, 0.25, 0.25],Q)

def K80(kappa = config['K80']['kappa']):
    ''' K80 '''
    return(TN93(alpha_1 = kappa, alpha_2 = kappa, beta = 1))


def LG():
    '''LG'''
    return(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],\
    [0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.0646, 0.022951, 0.042302, 0.04404, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147],\
    [[-16.060388, 0.425093, 0.276818, 0.395144, 2.489084, 0.969894, 1.038545, 2.06604, 0.358858, 0.14983, 0.395337, 0.536518, 1.124035, 0.253701, 1.177651, 4.727182, 2.139501, 0.180717, 0.218959, 2.54787], [0.425093, -21.193959, 0.751878, 0.123954, 0.534551, 2.807908, 0.36397, 0.390192, 2.426601, 0.126991, 0.301848, 6.326067, 0.484133, 0.052722, 0.332533, 0.858151, 0.578987, 0.593607, 0.31444, 0.170887], [0.276818, 0.751878, -17.692284, 5.076149, 0.528768, 1.695752, 0.541712, 1.437645, 4.509238, 0.191503, 0.068427, 2.145078, 0.371004, 0.089525, 0.161787, 4.008358, 2.000679, 0.045376, 0.612025, 0.083688], [0.395144, 0.123954, 5.076149, -24.184967999999994, 0.062556, 0.523386, 5.24387, 0.844926, 0.927114, 0.01069, 0.015076, 0.282959, 0.025548, 0.017416, 0.394456, 1.240275, 0.42586, 0.02989, 0.135107, 0.037967], [2.489084, 0.534551, 0.528768, 0.062556, -13.058626, 0.084808, 0.003499, 0.569265, 0.640543, 0.320627, 0.594007, 0.013266, 0.89368, 1.105251, 0.075382, 2.784478, 1.14348, 0.670128, 1.165532, 1.959291], [0.969894, 2.807908, 1.695752, 0.523386, 0.084808, -18.724285999999996, 4.128591, 0.267959, 4.813505, 0.072854, 0.582457, 3.234294, 1.672569, 0.035855, 0.624294, 1.223828, 1.080136, 0.236199, 0.257336, 0.210332], [1.038545, 0.36397, 0.541712, 5.24387, 0.003499, 4.128591, -23.8339, 0.348847, 0.423881, 0.044265, 0.069673, 1.807177, 0.173735, 0.018811, 0.419409, 0.611973, 0.604545, 0.077852, 0.120037, 0.245034], [2.06604, 0.390192, 1.437645, 0.844926, 0.569265, 0.267959, 0.348847, -16.613733, 0.311484, 0.008705, 0.044261, 0.296636, 0.139538, 0.089586, 0.196961, 1.73999, 0.129836, 0.268491, 0.054679, 0.076701], [0.358858, 2.426601, 4.509238, 0.927114, 0.640543, 4.813505, 0.423881, 0.311484, -6.179003, 0.108882, 0.366317, 0.697264, 0.442472, 0.682139, 0.508851, 0.990012, 0.584262, 0.597054, 5.306834, 0.119013], [0.14983, 0.126991, 0.191503, 0.01069, 0.320627, 0.072854, 0.044265, 0.008705, 0.108882, -22.217359, 4.145067, 0.159069, 4.273607, 1.112727, 0.078281, 0.064105, 1.033739, 0.11166, 0.232523, 10.649107], [0.395337, 0.301848, 0.068427, 0.015076, 0.594007, 0.582457, 0.069673, 0.044261, 0.366317, 4.145067, -31.146227, 0.1375, 6.312358, 2.592692, 0.24906, 0.182287, 0.302936, 0.619632, 0.299648, 1.702745], [0.536518, 6.326067, 2.145078, 0.282959, 0.013266, 3.234294, 1.807177, 0.296636, 0.697264, 0.159069, 0.1375, -19.90551, 0.656604, 0.023918, 0.390322, 0.748683, 1.136863, 0.049906, 0.131932, 0.185202], [1.124035, 0.484133, 0.371004, 0.025548, 0.89368, 1.672569, 0.173735, 0.139538, 0.442472, 4.273607, 6.312358, 0.656604, -24.724051, 1.798853, 0.099849, 0.34696, 2.020366, 0.696175, 0.481306, 1.898718], [0.253701, 0.052722, 0.089525, 0.017416, 1.105251, 0.035855, 0.018811, 0.089586, 0.682139, 1.112727, 2.592692, 0.023918, 1.798853, -23.599417, 0.094464, 0.361819, 0.165001, 2.457121, 7.803902, 0.654683], [1.177651, 0.332533, 0.161787, 0.394456, 0.075382, 0.624294, 0.419409, 0.196961, 0.508851, 0.078281, 0.24906, 0.390322, 0.099849, 0.094464, -16.74927, 1.338132, 0.571468, 0.095131, 0.089613, 0.296501], [4.727182, 0.858151, 4.008358, 1.240275, 2.784478, 1.223828, 0.611973, 1.73999, 0.990012, 0.064105, 0.182287, 0.748683, 0.34696, 0.361819, 1.338132, -12.953694999999998, 6.472279, 0.248862, 0.400547, 0.098369], [2.139501, 0.578987, 2.000679, 0.42586, 1.14348, 1.080136, 0.604545, 0.129836, 0.584262, 1.033739, 0.302936, 1.136863, 2.020366, 0.165001, 0.571468, 6.472279, -20.294897, 0.140825, 0.245841, 2.188158], [0.180717, 0.593607, 0.045376, 0.02989, 0.670128, 0.236199, 0.077852, 0.268491, 0.597054, 0.11166, 0.619632, 0.049906, 0.696175, 2.457121, 0.095131, 0.248862, 0.140825, -19.666723, 3.151815, 0.18951], [0.218959, 0.31444, 0.612025, 0.135107, 1.165532, 0.257336, 0.120037, 0.054679, 5.306834, 0.232523, 0.299648, 0.131932, 0.481306, 7.803902, 0.089613, 0.400547, 0.245841, 3.151815, -20.084989000000004, 0.249313], [2.54787, 0.170887, 0.083688, 0.037967, 1.959291, 0.210332, 0.245034, 0.076701, 0.119013, 10.649107, 1.702745, 0.185202, 1.898718, 0.654683, 0.296501, 0.098369, 2.188158, 0.18951, 0.249313, -18.608258]]\
    )

def loglikelihood_substitution(t,Q,eqFreq,alphabet,tandem_repeat):

    """ Calculate the likelihood of a repeat assuming a star tree and a sequence model
        defined by ``Q``, ``t``, ``eqFreq`` and ``alphabet``.

    Args:
        t (float): The divergence of the ``tandem_repeat`` units.
        Q (dict ?): A substitution matrix.
        eqFreq (dict ?): Equilibrium frequencies of all chars.
        alphabet (dict of str,?):  An alphabet from a substitution matrix.
        tandem_repeat (Repeat): An instance of the ``Repeat`` class.

    Returns:
        float: The loglikelihood value for the tandem repeat.

    """

    P = scipy.linalg.expm(Q*t)

    # msaTDN is a list of tandem_repeat.n numpy arrays.
    # Each numpy array represents the elements of the alphabet. The entries are the count of an alphabet element in the respective tandem_repeat column
    try:
        tandem_repeat.msaTDN
    except:
        tandem_repeat.msaTDN = []
        for column in tandem_repeat.msaTD_standard_aa:
            column = column.replace("-","")
            my_diag = {i:0 for i in range(len(alphabet.keys()))}
            for iN in column:
                my_diag[alphabet[iN]] = my_diag[alphabet[iN]] + 1
            tandem_repeat.msaTDN.append(np.array([j for j in my_diag.values()]))

    # First sum: Over all sites
    # Second sum: Over all possible ancestral characters
    # Third product: Over all repeat units
    likelihood = sum( np.log( np.sum( iJ * np.product(np.power(P[i], column)) for i,iJ in enumerate(eqFreq) ) ) for column in  tandem_repeat.msaTDN )
    return likelihood

def loglikelihood_gaps_starphylogeny_zipfian(t, tandem_repeat,
        indelRatePerSite = config['indel'].as_float('indelRatePerSite'),
        gaps = config['indel']['gaps'], indelZipf =config['indel'].as_float('zipf')):

    """ Calculate the log-likelihood of the gap structure in a ``tandem repeat`` assuming
     a star tree.

    Calculate the log-likelihood of the gap structure in a ``tandem repeat`` assuming a
    star tree (i.e. repeat unit correlations are not modeled). The gap length model is a
    zipfian distribution with parameter ``indelZipf``. The gap generation model is a
    simple poisson model with rate parameter ``indelRatePerSite``.

    Args:
        t (float): The divergence of the ``tandem_repeat`` units.
        tandem_repeat (Repeat): An instance of the ``Repeat`` class.
        indelRatePerSite (float): The mutation rate of insertions an deletions, as opposed
        to ``t``.
        gaps (str): The mode of gap counting. Options are "row_wise",
        "ignore_coherent_deletions", "ignore_trailing_gaps" and
        "ignore_trailing_gaps_and_coherent_deletions".
        indelZipf (float): The parameter of the Zipfian distribution.

    Returns:
        float: The loglikelihood value for the gap structure of the tandem repeat.

    ..  todo:: USE LOG EARLIER, NOT JUST IN RETURN STATEMENT
    """

    insertions = tandem_repeat.insertions
    deletions = tandem_repeat.deletions[gaps]
    gaps = tandem_repeat.gaps[gaps]

    # Calculate likelihood of the gap lengths (insertions and deletions combined)
    # Here, we use the Zipfian distributions
    l_length = scipy.special.zeta(indelZipf,1) ** - len(gaps)
    for gap in gaps:
        l_length = l_length * (gap ** -indelZipf)

    # Calculate the likelihood for the ratio of columns were an insertion has taken place,
    # as opposed to the number of columns, where no insertion has taken place.
    # Do the same for deletions.
    # Here, we use the Poisson model
    probability_gap_per_site = t * indelRatePerSite / 2
    l_insertions = scipy.stats.distributions.binom.pmf(len(insertions), tandem_repeat.lD * tandem_repeat.n + 1, probability_gap_per_site)
    l_deletions = scipy.stats.distributions.binom.pmf(len(deletions), tandem_repeat.lD * tandem_repeat.n, probability_gap_per_site)

    # Return the total likelihood of tandem_repeat's gap structure:
    return  np.log(l_length * l_insertions * l_deletions)

def loglikelihood_substitutions_gaps(t, substitution_args, gap_args):
    return loglikelihood_substitution(t, *substitution_args) + \
        loglikelihood_gaps_starphylogeny_zipfian(t, *gap_args)


def loglikelihood_starTopology_local(tandem_repeat,
            evolution_model=config['evolutionary_model'],
            gaps = config['indel']['gaps'],
            indelRatePerSite = config['indel'].as_float('indelRatePerSite'),
            parameters = False):

    if tandem_repeat.lD == 0:
        return 100,-100

    if parameters:
        Q,eqFreq,alphabet = parameters
    else:
        Q,eqFreq,alphabet = loadModel(evolution_model)

    if gaps:
        divergence, loglikelihoodDuplicationHistory = optimisation(function = loglikelihood_substitutions_gaps, args = [[Q,eqFreq,alphabet,tandem_repeat],[tandem_repeat,indelRatePerSite, gaps]])
    else:
        divergence, loglikelihoodDuplicationHistory = optimisation(function = loglikelihood_substitution, args = [Q,eqFreq,alphabet,tandem_repeat])

    return divergence,loglikelihoodDuplicationHistory

def loglikelihood_random(repeat,
                        evolution_model = config['evolutionary_model'],
                        parameters = False):

    if repeat.sequence_type == 'AA':
        if parameters:
            equilibrium_freq = parameters
        else:
            aaRatefilename = join(DATAROOT, "substitution_rate_matrices", evolution_model + ".dat") # load equilibrium frequencies
            equilibrium_freq = load_equilibrium_freq(aaRatefilename)

        loglikelihoodRandom = 0.
        for aa, freq in equilibrium_freq.items():
            # np.log uses the natural logarithm.
            loglikelihoodRandom += repeat.textD_standard_aa.count(aa) * np.log(float(freq))
    elif repeat.sequence_type == 'DNA':
        ## CHECK THIS!
        ## At the moment, we assume equal frequencies for each nucleotide - Baseml Model J69 or K80
        # python2: float for division neccesary  # np.log uses the natural logarithm.
        loglikelihoodRandom = repeat.totD * np.log(0.25)
    return loglikelihoodRandom

def phyloStarTopology_local(tandem_repeat, evolution_model=False, gaps = False, indelRatePerSite = 0.001, parameters = False):
    if not evolution_model:
        if tandem_repeat.sequence_type == 'AA':
            evolution_model = 'lg'
        else:
            evolution_model = 'k80'

    divergence, logLH1 = loglikelihood_starTopology_local(tandem_repeat=tandem_repeat, evolution_model=evolution_model, gaps = gaps, indelRatePerSite = indelRatePerSite, parameters = parameters)
    logLH0  = loglikelihood_random(repeat=tandem_repeat, evolution_model=evolution_model)

    return divergence, 2*(logLH1-logLH0)


###################### TANDEM REPEAT SCORE CALIBRATION ###################################

def calc_score(repeats, result_path, result_filename,
            scoreslist = config['score_calibration'].as_list('scoreslist'),
            save_calibration = config['score_calibration'].as_bool('save_calibration'),
            precision = config['score_calibration'].as_int('precision') ):

    sequence_type = repeats[0].sequence_type
    if repeats[0].sequence_type == 'AA':
        evolution_model = 'lg'
    else:
        evolution_model = 'k80'
    Q,eqFreq,alphabet = loadModel(evolution_model = evolution_model)
    parameters = [Q,eqFreq,alphabet]

    if sequence_type == 'AA':
        aaRatefilename = join(DATAROOT, "substitution_rate_matrices", evolution_model + ".dat") # load equilibrium frequencies
        equilibrium_freq = load_equilibrium_freq(aaRatefilename)
    else:
        equilibrium_freq = {i: 0.25 for i in ['A','C','G','T']}

    random = [loglikelihood_random(repeat=iRepeat, evolution_model=evolution_model, parameters = equilibrium_freq) for iRepeat in repeats]
    for iScore in scoreslist:
        if iScore == 'phylo':
            testStatistic = [loglikelihood_starTopology_local(iRepeat, gaps = False, parameters = parameters)[1] for iRepeat in repeats]
            testStatistic = [2*(iT-iR) if  iT != -1 else -1000 for iT,iR in zip(testStatistic,random)]
        elif iScore == 'phylo_gap01':
            testStatistic = [loglikelihood_starTopology_local(iRepeat, gaps = True, indelRatePerSite = 0.01, parameters = parameters)[1] for iRepeat in repeats]
            testStatistic = [2*(iT-iR) if  iT != -1 else -1000 for iT,iR in zip(testStatistic,random)]
        elif iScore == 'phylo_gap001':
            testStatistic = [loglikelihood_starTopology_local(iRepeat, gaps = True, indelRatePerSite = 0.001, parameters = parameters)[1] for iRepeat in repeats]
            testStatistic = [2*(iT-iR) if  iT != -1 else -1000 for iT,iR in zip(testStatistic,random)]
        elif iScore == 'phylo_paml':
            testStatistic = [loglikelihood_starTopology_paml(iRepeat) for iRepeat in repeats]
            testStatistic = [2*(iT-iR) if  iT != -1 else -1000 for iT,iR in zip(testStatistic,random)]
        elif iScore == 'parsimony':
            testStatistic = [np.round_(meanSimilarity(iRepeat, parsimony), decimals = precision) if iRepeat.lD != 0 else -1 for iRepeat in repeats]
        elif iScore == 'pSim':
            testStatistic = [np.round_(meanSimilarity(iRepeat, pSim), decimals = precision) if iRepeat.lD != 0 else -1 for iRepeat in repeats]
        elif iScore == 'entropy':
            testStatistic = [meanSimilarity(iRepeat, entropy) if iRepeat.lD != 0 else -1 for iRepeat in repeats]

        if save_calibration:
            filename = os.path.join(result_path,iScore,result_filename)
            print(filename)
            np.savez(filename , np.sort(testStatistic))
        else:
            if not hasattr(repeats[0],'score'):
                for iR in repeats:
                    iR.score = {}
            for iR,iT in zip(repeats,testStatistic):
                iR.dScore[iScore] = iT

def calibrate_score(repeats, resultFilePath,
                    scoreslist=config['score_calibration'].as_list('scoreslist')):

    for iScore in scoreslist:
        testStatistic = np.sort([iRepeat.dScore[iScore] for iRepeat in repeats])
        np.savez(os.path.join('_'.join([resultFilePath,iScore])) , testStatistic)

#############  CALCULATE AND SAVE THE PDF OF SCORES (Potentially deprecated) #############
def save_distribution(values, items, resultFilePath, fileName, inverse):
    for iC in range(len(items)):
        val = np.array([i[iC] for i in values])
        result_val, inv_ndx = np.unique(val, return_inverse=True)
        result_p = np.bincount(inv_ndx)
        result_p = result_p/float(sum(result_p))
        if inverse and items[iC] in ['phylo', 'pSim']: ## Are high values the best?
            result_p = result_p[::-1]
            result_val = result_val[::-1]
        result_p = result_p.cumsum()
        if not os.path.isdir(os.path.join(resultFilePath, items[iC])):
            os.makedirs(os.path.join(resultFilePath, items[iC]))
        np.savez(os.path.join(resultFilePath, items[iC], fileName), cdf=result_p,value=result_val)

def calculatePDFScores(repeats, resultFilePath, fileName,
                scoreslist = config['score_calibration'].as_list('scoreslist')):

    """ CALCULATE THE PROBABILITY DISTRIBUTION OF SCORES """

    # Can you generalise the next command for arbitrary classifiers?
    scores = [[repeat.score(scoreslist[0]), repeat.score(scoreslist[1])] for repeat in repeats]
    save_distribution(values = scores, items = scoreslist, resultFilePath = resultFilePath, fileName = fileName, inverse = True)
