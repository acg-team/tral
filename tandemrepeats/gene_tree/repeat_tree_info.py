# (C) 2013 Elke Schaper

import itertools,logging,os,random,re,sys
import numpy as np
import scipy as sc
from scipy.misc import comb
from scipy.stats import hypergeom

from repeat.paths import *

sys.path.append(os.path.join(ROOT,'lib/python3.2.2/'))
from ete2 import Tree

# Michiel de Hoon's PyCluster: http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm
import Bio.Cluster

logger = logging.getLogger(__name__)


################################ Repeat Tree class #######################################


class Repeat_Tree(Tree):
            

    def subtree(self, leaves):
        
        self.unroot()
        # Remove leaves that not in <leaves>
        leaves_remove = [iLeaf for iLeaf in self.iter_leaves() if not iLeaf.name in leaves]
        self.remove_leaves(leaves_remove)

    def remove_leaves(self,leaves):
    
        ''' remove <leaves> and afterwards unnecessary parent nodes from <tree>
            Warning: What happens if the parent of a node is the tree's root?
        '''
        
        for iLeaf in leaves:
            parent = iLeaf.up
            if parent.is_root():
                # Define a different root, as the current will be deleted
                # Choose a closest sister as the ungroup. 
                sister = iLeaf.get_sisters()[0]
                self.set_outgroup(sister)
                self.unroot()
                parent = iLeaf.up
        
            iLeaf.delete()
            parent.delete()

    def set_Ensembl_ID_and_index(self):
        lID = []
        for iLeaf in self.iter_leaves():
            ensembl_ID, index = iLeaf.name.split("_")
            iLeaf.add_features(ensembl_ID=ensembl_ID)
            iLeaf.add_features(index=index)
            lID.append(ensembl_ID)
        self.lID = set(lID)
        
        # Calculate the number of repeat units for each ensembl ID:
        self.n = {i: lID.count(i) for i in self.lID}

    def get_cherries(self):
    
        lLeaf = list(self.iter_leaves())
        
        if len(self) <= 3:
           lCherry = list(itertools.permutations(lLeaf))
        
        self.lCherry = []
        while lLeaf:
            iL = lLeaf.pop()
            lSister = iL.get_sisters()
            for iS in lSister:
                if iS.is_leaf():
                    self.lCherry.append((iL, iS))
                    lLeaf.remove(iS)
        
        self.nC = len(self.lCherry)
                
    def get_bisample_cherries(self):
    
        if not hasattr(self, 'lCherry') or not hasattr(self, 'nC'):
            self.get_cherries()
        if not hasattr(self, 'n'):
            self.set_Ensembl_ID_and_index()
            
        self.lBisample_cherry = [iCH for iCH in self.lCherry if iCH[0].ensembl_ID != iCH[1].ensembl_ID]
        
        # Save the number of bisample cherries
        self.nBC = len(self.lBisample_cherry)
        
        # save the p-Value of <self.nBC> given <self.nC>, and <self.n>
        n = list(self.n.values())
        if len(n) != 2:
            logging.debug("There are {} ensembl identifiers in this repeat units tree [Expect: 2].".format(len(n)))
        
        if not hasattr(self, 'pValue'):
            self.pValue = {}
        
        if n[0] <= 1 or n[1] <= 1:
            self.pValue['nBC'] = 1  
        else:
            self.pValue['nBC'] = p_Value_bisample_cherries(n = n,n_cb = self.nBC,n_c = self.nC)
        
        
    def get_umbel_pairs(self):
    
        '''  Collapse all umbels to one leaf. Then, return cherries 
            To start, define the root in a smart way (e.g. make sure the root is not in
            the middle of an umbel.)
        '''
        
    
    def get_min_split(self):
        
        ''' What is the minimal number of splits to separate repeat units from both TRs? '''
        
        cost = get_parsimony_cost_state(self,list(self.lID))
        self.min_split = min(cost.values())
            
    
    def order_conservation_test(self, sample_type = 'bisample_cherries', test = 'Kendall', nPermutation = 1000):
    
        if sample_type == 'bisample_cherries':
            sample = self.lBisample_cherry
        
        if test == 'Spearman':
            test_function = spearman
        elif test == 'Kendall':
            test_function = kendall
        
        if not hasattr(self, 'order_conversation'):
            self.order_conservation = {}
        
        if len(sample) < 2:
            self.order_conservation[test] = None
        else:
            # Extract the repeat unit index.
            index_pairs = [[iS[0].index,iS[1].index] for iS in sample]  
        
            # Transpose the repeat unit indices.
            index_pairs_t = list(map(list, zip(*index_pairs)))
    
            test_statistic = test_function(index_pairs_t)
            shuffle_test = []
            for count in range(nPermutation):
                random.shuffle(index_pairs_t[0])
                shuffle_test.append(test_function(index_pairs_t))
         
            pValue = len([i for i in shuffle_test if i >= test_statistic])/nPermutation

            self.order_conservation[test] = {'test_statistic':test_statistic, 'pValue':pValue, 'nPermutation':nPermutation}
            
################################# ORDER CONSERVATION #####################################

def spearman(index_pairs_t):
    return 1 - Bio.Cluster.distancematrix((index_pairs_t[0],index_pairs_t[1]), dist="s")[1][0]
    
def kendall(index_pairs_t):
    return 1 - Bio.Cluster.distancematrix((index_pairs_t[0],index_pairs_t[1]), dist="k")[1][0]

################################# Minimal number of splits ###############################


def get_parsimony_cost_state(tree_node, lState):

    if tree_node.is_leaf():
        return {iState: (0 if tree_node.ensembl_ID == iState else 1) for iState in lState}
        
    else:
        cost_children = [get_parsimony_cost_state(iC, lState) for iC in tree_node.get_children()]
        my_cost = {iState: sum(iC[iState] for iC in cost_children) for iState in lState}
        my_cost[lState[0]] = min(my_cost[lState[0]],my_cost[lState[1]] + 1)
        my_cost[lState[1]] = min(my_cost[lState[1]],my_cost[lState[0]] + 1)
        return my_cost
     
############################### Number of bisample cherries ##############################
        
             
def p_Value_bisample_cherries(n,n_cb,n_c):
    
    ''' Return the probability of <n_bc> or more bisample cherries with <n> tips of either
        colour and <n_c> cherries in total.
        
        Parameters e.g.:
        n = [5,9]
        n_bc = 3
        n_c = 4
        would return ~0.6873
        
        At current, do not check whether inputs are sensible.
    
    '''
    
    p = 0
    for i in range(n_cb, min(n[0],n[1],n_c) + 1):
        p += p_bisample_cherries(n,i,n_c)
    return p


#def p_Value_bisample_cherries_2(n,n_cb,n_c):
#    
#    return comb(n_c,n_cb) * np.power(2,n_cb) * np.prod(range(n[0]-n_cb+1, n[0]+1)) * np.prod(range(n[1]-n_cb+1, n[1]+1)) / np.prod(range(n[0]+n[1]-2*n_cb+1, n[0]+n[1] + 1))


def n_0_bisample_cherries(n,n_c):
    
    '''
    Return the # of possibilities of exactly 0 bisample cherries with <n> tips of either
    colour and <n_c> cherries in total.
    '''   
    
    count = 0
    for n_c0 in range( max(0, int(n_c - np.floor(n[1]/2))) , min(n_c, int(np.floor(n[0]/2))) +1 ):
        count += comb(n_c,n_c0) * comb(np.sum(n) - 2*n_c, n[0] - 2*n_c0)
    return count


def p_bisample_cherries(n,n_cb,n_c):
    
    '''
    Return the probability of exactly <n_bc> bisample cherries with <n> tips of either
    colour and <n_c> cherries in total.
    
    Exp test to check correctness:
    ###
    import random
    order = [0]*n[0] + [1]*n[1]
    nTest = 10000
    exp = {i:0 for i in range(n_c+1)}
    for i in range(nTest):
        random.shuffle(order)
        exp[sum([0 if order[2*j] == order[2*j+1] else 1 for j in range(n_c)])] += 1
    exp = {i:j/nTest for i,j in exp.items()}
    ###
    
    '''
    return  np.power(2,n_cb) * n_0_bisample_cherries([i-n_cb for i in n], n_c - n_cb) * comb(n_c,n_cb) / comb(np.sum(n),n[0])


########## PROBABILITY OF EXTREME CASE OF REPEAT UNIT PHYLOGENY ON RANDOM TREE ###########

#  GET RANDOM n0 n1 DATA
# import numpy as np
# import scipy as sc
# from scipy.stats.distributions import binom as binom
# import random
# from collections import defaultdict
# # Input: Available n_0,n_1 combinations
# p_combinations = defaultdict(int)
# n_min = 2
# n_max = 5
# nTests = 1000
# for i in range(nTests):
#     tmp_int = random.randint(n_min,n_max)
#     p_combinations[(tmp_int,random.randint(tmp_int,tmp_int+1))] += 1
# 
# 
# # Take only combinations where n0 = n1
# p_possibly_perfect = {i[0]:j for i,j in p_combinations.items() if i[0]==i[1]}
# 
# 
# # Get probability for each value
# p_distributions = [[binom.pmf(iCount,count,probability_perfect(n)) for iCount in range(count)] for n,count in p_possibly_perfect.items()]
# 
# # Calculate convolved distribution
# p_final = p_distributions[0]
# for i in p_distributions[1:]:
#     p_final = sc.convolve(p_final, i)
# 
# # Calculate inversed cumulative distribution
# p_final = sc.cumsum(p_final[::-1])
# 
# # ES: CHECK FINAL STEPS.
# # Calculate p-Value 
# k = 100
# p_final[k]

def probability_perfect_conservation(n):
        
    ''' Return the probability of pairwise perfect TR unit conversation occuring by chance.
        Assume that phylogenies are uniformly distributed.
        
        Input: Number of TR units in each TR n
        I.e. The total number of TR units is 2*n. '''
    if n >= 50:
        print("Warning! The resulting probability is smaller than the precision of the system for n = {}".format(n))
    return np.power(2,n) /sc.misc.factorial(4*n-4) * sc.misc.factorial(2*n-2) * sc.misc.factorial(2*n-4)/sc.misc.factorial(n-2)
        
def probability_completely_diverged(n0, n1):
    ''' Return the probability of complete TR unit separation (parsimony value / min_split = 1) occuring by chance.
        Assume that phylogenies are uniformly distributed.
        
        Input: Number of TR units in each TR n0 and n1
        I.e. The total number of TR units is n0 + n1. '''
    return fSteel(n0, n1, 1)

def fSteel(a, b, k):
      
    return sc.misc.factorial(k - 1) * (2*(a + b) - 3*k) /bSteel(a + b - k + 2)  * NSteel(a, k) * NSteel(b, k)

def NSteel(n, k):
    
    if n < k:
        return 0
    else:
        return sc.misc.factorial(2*n - k - 1)/ sc.misc.factorial(n - k) / sc.misc.factorial(k - 1) / np.power(2, n - k)

def bSteel(n):
    return sc.misc.factorial2(2*n - 5)

################################ DUPLICATION HISTORY (probably outdated) #################

def detect_speciation_duplication_events(tree):
    # Traverse phylogenetic tree to mark all speciation ('S') and duplication ('D') nodes
    unvisited_nodes = [i for i in tree.iter_leaves()]
    while len(unvisited_nodes)>0:
        brother = unvisited_nodes.pop(0)
        sister = brother.get_sisters()[0]
        if hasattr(sister,'ensembl_ID'):
            try:
                unvisited_nodes.remove(sister)
            except:
                pass
            parent = brother.up
            parent.add_features(ensembl_ID = brother.ensembl_ID | sister.ensembl_ID)
            if not parent.is_root():
                unvisited_nodes.append(parent)
            if parent.ensembl_ID == brother.ensembl_ID or parent.ensembl_ID == sister.ensembl_ID:
                # Duplication event.
                parent.add_features(type = 'D')
            else:
                # Speciation event.
                parent.add_features(type = 'S')
        else:
            unvisited_nodes.append(brother)
       
    ## Get speciation nodes:
    speciation_nodes = [iNode for iNode in tree.iter_descendants() if iNode.type == 'S']   
    return speciation_nodes

def duplication_history_test():

    ''' IDEA: Make a dictionary based check to see what duplication history events explain the repeat unit order in a phylogeny. 
    NOT COMPLETED. '''
    ancestor = tree.get_common_ancestor("ENSONIP00000010035_4", "ENSONIP00000010035_9")
    
    # Reorder children of ancestor 
    # Child with the greatest #leaves is considered the first child.
    # (Later: When #leaves is equal, topology information should be included into the ordering.)
    
    for iNode in ancestor.traverse():
        # Inverse the order of the child nodes if the second node has more leaves than the first.
        if len(iNode.children) == 2 and len(iNode.children[0]) < len(iNode.children[1]):
            iNode.children = [iNode.children[1],iNode.children[0]]
        ###  CHECK! ete2 ladderize might do the same or something better!
    
    
    # Get the tree topology in Newick format (lacking all name and distance information)
    sub_tree_topology = ancestor.write(format=100)[:-1]
    # Get the leaves' names in the same order
    sub_tree_leaves = [i.name for i in ancestor.iter_leaves()]
    # Translate the leaf names to repeat unit numbers.
    sub_tree_coordinates = [leaves[i][1] for i in sub_tree_leaves]
    
    
    # Collapse neighbouring nodes
    
    
    duplication_topologies[sub_tree_topology]
    
    duplication_topologies = {'((,),)': {'D': {1:1,2:1} , 'R':{1:1}},
                              '((,),(,))': {'D': {1:1,2:1}},
                              '((,),),)': {'1324': {'D': {1:2,2:1},'R':{1:1}}, '1342':{'D': {1:2,3:1},'R':{1:2}}, '1423':{'D': {1:2,3:1},'R':{2:1}}},
                              '((,),(,),)': {'13245':{'D':{1:2,2:1}}, '13254':{'D':{1:2,2:1,5:1},'R':{1:1,1:3}}, '14253':{'D':{1:2,3:1},'R':{1:1}}, '15243':{'D':{1:2,2:1,5:1},'R':{1:1,2:2}}}, #alternative for the last entry: 'R':{1:2,3:1}
                            }

