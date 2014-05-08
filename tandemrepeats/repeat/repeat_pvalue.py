# (C) 2012, Elke Schaper

import os
import bisect
from os.path import join
import csv
import logging
import numpy as np
import scipy as sp
import scipy.stats, scipy.special

logger = logging.getLogger(__name__)

from repeat.paths import *

path_score = join(DATAROOT, 'pvalue')

########################## REPEAT SCORE p-VALUE CALCULATION FUNCTIONS ####################

##################### phylo & entropy ############################

def empiricalList(l,n, sequence_type = 'AA', score = 'phylo'):

    ''' loads and returns a numpy list with 10,000 empirical values
        of the user defined distribution from memory. 
        If no values are known for l and n, the closest values for l and n are chosen instead.'''
        
    lMax = 99
    nMax = 49
    if score == 'entropy':
        lMax = 20
        nMax = 5

    totalRepeatLengthMax = 1000
    l = min(l, lMax)
    n = min(n, nMax)
    if l*n > totalRepeatLengthMax:
        print("l: %d and n: %d are bigger than the totalRepeatLengthMax: %d" %(l,n,totalRepeatLengthMax))
        ## Dirty little hack: as we do not have data for this pair of l and n, apply nearest data file.
        while l*n > totalRepeatLengthMax:
            if l>n:
                l -= 1
            else:
                n -= 1
    
    file = join(path_score, sequence_type, score, str(l) + '_' + str(n) + '.npz')
    if not os.path.isfile(file):
        raise ValueError("complete pdf file %s does not exist!" % file)  
    
    myEmpiricalList = np.load(file)
    ## It is absolutely necessary to close the numpy filehandle.
    ## Otherwise, you will have too many open operating system filehandles
    ## if you run this function many times (more than ulimit -n allows, that is)
    empiricalList = myEmpiricalList['arr_0']
    myEmpiricalList.close()    
    
    return empiricalList

def pValueFromEmpiricialList(myTR, score = 'phylo', myScore = None, empirical = []):
    
    '''  Calculates the p-Value of a score for the given myTR. The p-Value is the number of 
    scores for comparable tandem repeats, i.e. of same repeat unit length and repeat unit copy number
    that are as good or better. 
    '''    
    
    if myScore == None:
        myScore = myTR.score[score]
        
    if len(empirical) == 0:
        empirical = empiricalList(myTR.lD, myTR.n, myTR.sequence_type, score)   
    
    if score in ['entropy']: # A smaller score is a better score for all scores in this list.
        if myScore == -1:
            return 1
        else:
            return bisect.bisect_right(empirical, myScore)/len(empirical)
    else:
        return 1 - bisect.bisect_left(empirical, myScore)/len(empirical)
  
    

########################### pSim & parsimony: read in PDF ###########################

def columnPDF(n, score = 'psim', sequence_type = 'AA'):
    ''' Digs aus the pdf of the score on random sequence_type data of length n
        The order of the pdf is WORST FIRST!'''

    # Currently, the pdfs are only save up to nMax = 150
    nMax = 150
    with open(join(path_score, sequence_type, score, str(min(n,nMax)) + '.txt'), 'r') as pdf_file:
        pdf = [i for i in csv.reader(pdf_file, dialect='excel-tab')]
        pdf = [np.double(i[0]) for i in pdf[1:]]
        return np.array(pdf)

def dAverageMultinom(l,n, sequence_type, score):
    ## Returns the l-times self-convoluted pdf of the score on random sequence_type data of length n
    ## The order of the pdf is kept from columnPDF()

    ## score in ['psim', 'parsimony']
    ## sequence_type in ['AA', 'DNA']
    completePDF = columnPDF(n = n, score = score, sequence_type = sequence_type)
    if l != 1:
        singleColumnPDF = completePDF
        for i in range(2,l+1):
            completePDF = sp.convolve(completePDF, singleColumnPDF)
    return completePDF

######################### pSim ###########################

def calculate_repeat_structure(myTR):
    # You can use a different null distribution for each column of different length for both
    # the parsimony and the pSim score.
    # Therefore, calculate the structure of the TR before you calculate the null distribution
    # of tandem repeat scores of tandem repeats with the same gap distribution.
    
    repeatStructure = [
            len(column.replace("-",""))
            for column in myTR.msaTD if (2*column.count('-') < len(column))
    ]
    nRepeatStructure = list(set(repeatStructure))
    lRepeatStructure = [
        repeatStructure.count(i) for i in nRepeatStructure
    ]
    return lRepeatStructure,nRepeatStructure
    


def dAverageMultipleMaxMultinom(myTR, precision = 10000.): # python 2: precision must be float
    # if precision higher than max(uint32) use uint64 instead
    # CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html
    

    lRepeatStructure,nRepeatStructure = calculate_repeat_structure(myTR)
    
    p = dAverageMultinom(
        lRepeatStructure[0], nRepeatStructure[0],
        myTR.sequence_type, score = 'psim'
    )
    val = np.array(
        (np.r_[lRepeatStructure[0]:(lRepeatStructure[0] * nRepeatStructure[0] + 1.)]
         * (precision/nRepeatStructure[0])).round(),
        dtype='uint32'
    )
    if len(lRepeatStructure) == 1:
        # return best values first. for pSim, best = 1, worst = 0
        return  p[::-1], val[::-1]/(precision*sum(lRepeatStructure))
    else:
        for i in range(1,len(lRepeatStructure)):
            p = np.outer(p,
                dAverageMultinom(
                    lRepeatStructure[i], nRepeatStructure[i],
                    myTR.sequence_type, score = 'psim')
            ).ravel()

            val = (val[:,np.newaxis] +
                np.array(
                    (np.r_[lRepeatStructure[i]:(lRepeatStructure[i] * nRepeatStructure[i] + 1.)]
                        * (precision/nRepeatStructure[i])
                    ).round(),
                    dtype='uint32')
            ).ravel()


            x = np.bincount(val,weights=p)
            val = np.unique(val)
            p = x[val]
        # return best values first. for pSim, best = 1, worst = 0
        return p[::-1], (val/(precision*sum(lRepeatStructure)))[::-1]

def pValuePSim(myTR):
    precision = 10000.

    pdf = dAverageMultipleMaxMultinom(myTR,precision)
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(myTR.score['pSim'] == pdf[1])[0]
    index = np.where(np.abs(myTR.score['pSim'] - pdf[1]) <= 1./precision)[0]
    if len(index) == 1: ## standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(myTR.score['pSim'] > pdf[1])[0]
    if len(indices) >= 1:
    ## if pars is not exactly included in list, give back mean of value above
    ## and value below (should only occur due to numerical imprecision)
        a = min(indices)
        return((cumsumPDF[a] + cumsumPDF[a-1])/2)
    else: ## This should only occur if score['pSim'] is really bad or really good
        return 1

#################################### parsimony #################################


def dAverageMultipleParsMultinom(myTR, precision = 10000.): # python 2: precision must be float
    # if precision higher than max(uint32) use uint64 instead
    # CHECK: http://docs.scipy.org/doc/numpy/user/basics.types.html

    lRepeatStructure,nRepeatStructure = calculate_repeat_structure(myTR)    

    p = dAverageMultinom(
        lRepeatStructure[0], nRepeatStructure[0],
        myTR.sequence_type, score = 'parsimony'
    )
    val = np.array(
        (np.r_[(lRepeatStructure[0]*(nRepeatStructure[0]-1)):-1:-1.]
        * (precision/(nRepeatStructure[0]-1.))
        ).round(),
        dtype='uint32'
    )
    if len(lRepeatStructure) == 1:
        # return best values first. for parsimony, best = 0, worst = 1
        return  p[::-1], (val/(precision *sum(lRepeatStructure)))[::-1]
    else:
        for i in range(1,len(lRepeatStructure)):
            try:
                p = np.outer(p,
                    dAverageMultinom(
                        lRepeatStructure[i],
                        nRepeatStructure[i],
                        myTR.sequence_type, score = 'parsimony')
                ).ravel()

                val = (val[:,np.newaxis] + np.array(
                    (np.r_[(lRepeatStructure[i]*(nRepeatStructure[i]-1)):-1:-1.]
                     * (precision/(nRepeatStructure[i]-1.))
                    ).round(),
                    dtype='uint32')
                ).ravel()

                x = np.bincount(val,weights=p)
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
        return p, val/(precision*sum(lRepeatStructure))

def pValuePars(myTR):
    precision = 10000.

    pdf = dAverageMultipleParsMultinom(myTR,precision)
    if pdf == False:
        return 1
    cumsumPDF = np.cumsum(pdf[0])

    #index = np.where(pdf[1] == myTR.score['parsimony'])[0]
    index = np.where(np.abs(myTR.score['parsimony'] - pdf[1]) <= 1./precision)[0]
    if len(index) == 1: ## standard case: score is included in pdf[0]
        return(cumsumPDF[index[0]])

    indices = np.where(myTR.score['parsimony'] < pdf[1])[0]
    if len(indices) >= 1:
        ## if pars is not exactly included in list, give back mean of value
        ## above and value below (should only occur due to numerical imprecision)
        a = min(indices)
        return (cumsumPDF[a] + cumsumPDF[a-1])/2
    else: ## This should only occur if score['parsimony'] is really bad
        return 1



####################################### gap penalty #####################################

def gapPenalty(myTR, mu):
    return ((mu/(1-mu)) ** myTR.nGapStructure) * myTR.pGapStructure   
    
############################### pValue distribution ######################################

def calc_pValues(repeats, resultFilePath, fileName, scoreslist=['phylo','phylo_gap'], gappy_data = False):

    ''' save distribution of p-Values '''
    
    # If the repeats are not gappy, than you can use always the same distribution of scores
    # on random data to calculate the p-Value. Otherwise, there might be deletion columns
    # and the distribution should be loaded each time 'pValueFromEmpiricialList' is called.
    empirical_list = []

    for iScore in scoreslist:
        if 'parsimony' == iScore:
            testStatistic = [pValuePars(iRepeat) if iRepeat.lD != 0 else 1 for iRepeat in repeats]
        elif 'pSim' == iScore: 
            testStatistic = [pValuePSim(iRepeat) if iRepeat.lD != 0 else 1 for iRepeat in repeats]
        else:
            if not gappy_data:
                empirical_list = empiricalList(repeats[0].lD, repeats[0].n, repeats[0].sequence_type, iScore)
            testStatistic = [pValueFromEmpiricialList(iRepeat, iScore, empirical = empirical_list) if iRepeat.lD != 0 else 1 for iRepeat in repeats]
        score_path = os.path.join(resultFilePath,iScore)
        if not os.path.isdir(score_path):
                os.makedirs(score_path)
        print(os.path.join(score_path,fileName))        
        np.savez(os.path.join(score_path,fileName) , np.sort(testStatistic))        

