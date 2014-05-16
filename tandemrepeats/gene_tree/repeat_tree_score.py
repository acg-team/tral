import numpy as np
import scipy as sc
from scipy.stats.distributions import binom as binom
import random
from collections import defaultdict



# Input: Available n_0,n_1 combinations
p_combinations = defaultdict(int)
n_min = 2
n_max = 5
nTests = 1000
for i in range(nTests):
    tmp_int = random.randint(n_min,n_max)
    p_combinations[(tmp_int,random.randint(tmp_int,tmp_int+1))] += 1


# Take only combinations where n0 = n1
p_possibly_perfect = {i[0]:j for i,j in p_combinations.items() if i[0]==i[1]}


# Get probability for each value
p_distributions = [[binom.pmf(iCount,count,probability_perfect(n)) if n < 50 else [1]+[0]*(count-1) for iCount in range(count)] for n,count in p_possibly_perfect.items()]

# Calculate convolved distribution
p_final = p_distributions[0]
for i in p_distributions[1:]:
    p_final = sc.convolve(p_final, i)

# Calculate inversed cumulative distribution
p_final = sc.cumsum(p_final[::-1])

# ES: CHECK FINAL STEPS.
# Calculate p-Value 
k = 100
p_final[k]



def probability_perfect(n):
        
    ''' Return the probability of pairwise perfect TR unit conversation occuring by chance.
        Assume that phylogenies are uniformly distributed.
        
        Input: Number of TR units in each TR n
        I.e. The total number of TR units is 2*n. '''
    return np.power(2,n) /sc.misc.factorial(4*n-4) * sc.misc.factorial(2*n-2) * sc.misc.factorial(2*n-4)/sc.misc.factorial(n-2)
        
