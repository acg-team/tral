# (C) 2013 Elke Schaper
import numpy as np
import scipy as sc
from copy import copy
import logging, os, subprocess, tempfile
from Bio import AlignIO

logger = logging.getLogger('root')

from . import repeat_info
from repeat.paths import *

''' Implement a range of tests on the conversation of TR units between homologous TRs. 
The focus is on tests for the order of TR units. '''

################################## EDIT DISTANCES ########################################


def calc_pairwise_alignment_distance_matrix(msa_a, msa_b):

    ''' calculate lower left triangle matrix of distances pairwise alignment distances
        two lists of tandem repeat units A and B
        
        If you are interested in all pairwise distances, i.e. including the distances of the units within A or within B,
        provide msa_a+msa_b,msa_a+msa_b as parameter. The resulting distance matrix is quadratic.
        
        Parameters: MSAs for A and B.
    '''

    pairwise_alignment_distance = numpy.zeros(shape = (len(msa_a), len(msa_b)))
    
    for i,j in itertools.product(range(len(msa_a)),range(len(msa_b))):
        print('*')
        print(i)
        print(j)
        pairwise_alignment = ['','']
        # Remove gaps from the pairwise alignment
        for k,l in zip(msa_a[i], msa_b[j]):
            if (k != '-') or (l != '-'):
                pairwise_alignment[0] += k
                pairwise_alignment[1] += l
        # Calculate the distance for the current pair                
        pairwise_alignment = repeat_info.Repeat(pairwise_alignment)
        distance, score = repeat_score.phyloStarTopology_local(pairwise_alignment, gaps = 'row_wise', indelRatePerSite = 0.01)
        pairwise_alignment_distance[i][j] = distance
        
    return pairwise_alignment_distance


def create_random_pairwise_alignment_distance_matrix(n_a, n_b):
    
    n_all = n_a + n_b
    
    # Create n_allxn_all matrix. However, in practice, only lower left triangle of distance matrix is used.
    pairwise_alignment_distance = np.zeros(shape = (n_all,n_all))
    for i in range(1,n_all):
        for j in range(i):
            pairwise_alignment_distance[i][j] = np.random.rand()

    return pairwise_alignment_distance

def edit_distance_nj(pairwise_alignment_distance, n_a, n_b):

    ''' Calc the edit distance of the repeat units of two homologous tandem repeats A and B.
        
        Parameters: a <n_a>x<n_b> inverse distance matrix <pairwise_alignment_distance>
        The columns and rows of this matrix are ordered by
        first all repeat units is tr_a, and second, all repeat units in tr_b
        
        WARNING! This function expects the inverse distances as input.
        The highest value is treated as the lowest distance.
        
        Return the edit events'''
    
    active_nodes = list(range(n_a + n_b))
    # Initialise nodes 
    nodes = [{'tr_a': (i,i), 'tr_b': None} for i in range(n_a)] + [{'tr_a': None, 'tr_b': (i,i)} for i in range(n_b)]
    events = None
    
    while len(active_nodes) > 1:
        pairwise_alignment_distance, nodes, active_nodes, events = edit_distance_nj_iteration(pairwise_alignment_distance, nodes, active_nodes, events)
    
    return events


def edit_distance_nj_iteration(distance_matrix, nodes, active_nodes,events = None):
    
    # Initialise event count
    if events == None:
         events = {'duplication': 0, 'deletion': 0, 'translocation': 0, 'speciation': 0}
    
    # Find current max
    i_max = {'index': (0,0), 'value': 0}
    for index,j in enumerate(active_nodes):
        for i in active_nodes[:index]:
            if distance_matrix[j][i] > i_max['value']:
                i_max['value'] = distance_matrix[j][i]
                i_max['index'] = (j,i)
    
    # Update event count
    i = i_max['index'][1]
    j = i_max['index'][0]
    
    # It is i<j and leaves 0 come before leaves 1 in <distance_matrix>.
    # Check what type of event (speciation, duplication, translocation) occured.
    if nodes[i]['tr_b'] == None and nodes[j]['tr_a'] == None:
        # Speciation event
        events['speciation'] += 1
        # Update node[i] as merger of former node[i] and node[j]
        nodes[i] = {'tr_a': nodes[i]['tr_a'], 'tr_b': nodes[j]['tr_b']}
    elif (nodes[i]['tr_b'] == None and nodes[j]['tr_b'] == None) or (nodes[i]['tr_a'] == None and nodes[j]['tr_a'] == None):
        # Leaf duplication in 'tr_a' or 'tr_b'. Set <tr_d> to either.
        events['duplication'] += 1
        tr_d = 'tr_a' if nodes[i]['tr_b'] == None else 'tr_b'
        # Check for translocation:
        if not (nodes[i][tr_d][1] < nodes[j][tr_d][0] or nodes[i][tr_d][0] > nodes[j][tr_d][1]):
            events['translocation'] += 1
        # Update node[i] as merger of former node[i] and node[j]
        nodes[i][tr_d] = (min(nodes[i][tr_d][0],nodes[j][tr_d][0]), max(nodes[i][tr_d][1],nodes[j][tr_d][1]))
    else: #nodes[i]['tr_a'] == None or nodes[i]['tr_b'] == None or nodes[j]['tr_a'] == None or nodes[j]['tr_b'] == None: 
        # Duplication event 
        events['duplication'] += 1
        
        # Check for deletion:
        # An alternative to using <done> to break two loops at once would be to merge both for loops (e.g., itertools.product)
        done = False
        for node_d, node_n in zip([i,j], [j,i]):
            for tr_d,tr_n in zip(['tr_a','tr_b'],['tr_b','tr_a']):
                if nodes[node_d][tr_d] == None:
                    events['deletion'] += 1
                    # For the deletion tr <tr_d>, this merger contains the information from the non-deletion node <node_n> only.
                    nodes[i][tr_d] = nodes[node_n][tr_d]
                    #  Check for translocation in the none-deletion tr <tr_n>:
                    if not (nodes[i][tr_n][1] < nodes[j][tr_n][0] or nodes[i][tr_n][0] > nodes[j][tr_n][1]):
                        events['translocation'] += 1
                    # For the none-deletion tr <tr_n>, this merger contains the information from both nodes
                    nodes[i][tr_n] = (min(nodes[i][tr_n][0],nodes[j][tr_n][0]), max(nodes[i][tr_n][1],nodes[j][tr_n][1]))
                    done = True
                    break
            if done:
                break
        else:
            # No deletion has occurred:
            # Update node[i] as merger of former node[i] and node[j]
            for tr in ['tr_a','tr_b']:
                # Check for translocations
                if not (nodes[i][tr][1] < nodes[j][tr][0] or nodes[i][tr][0] > nodes[j][tr][1]):
                        events['translocation'] += 1
                nodes[i][tr] = (min(nodes[i][tr][0],nodes[j][tr][0]), max(nodes[i][tr][1],nodes[j][tr][1]))
    
    
    # Remove index from active nodes
    active_nodes.remove(j)    
        
        
    # Update  distance matrix        
    for k in active_nodes:
        if k < i:
            distance_matrix[k][i] = (distance_matrix[k][i]+distance_matrix[k][j])/2
        elif k == i:
            continue
        elif k > j:
            distance_matrix[i][k] = (distance_matrix[i][k]+distance_matrix[j][k])/2
        else:
            # i < k < j
            distance_matrix[i][k] = (distance_matrix[i][k]+distance_matrix[k][j])/2
    
    return distance_matrix, nodes, active_nodes, events


################################## AD HOC DISTANCES #####################################


def pairwise_order_conservation(a,b):

    ''' Count how often the order is conserved between pairwise best hits
     from two homologous tandem repeats A and B.
    
    Parameters: a: Index of best hit tandem repeat unit in B for each repeat unit in A
                b: vice versa.
                
    E.g. for perfect conservation:
    a = [0,1,2,3,4,5]
    b = [0,1,2,3,4,5]
    
    Direction: from A to B. I.e. For any pair of repeat units in A, can you predict the order
    of the corresponding best hit repeat units in B?
    
    Return the count of 
    - correct predictions ('conserved')
    - incorrect predictions ('not_conserved')
    - missing predictions in case of identical references from A to B ('identical')'''

    count = {'conserved': 0, 'identical': 0, 'not_conserved':0}

    for i in range(len(a)):
        for j in range(i):
            if a[i] == a[j]:
                count['identical'] += 1
            elif a[i] > a[j]:
                count['conserved'] += 1
            else:
                count['not_conserved'] += 1

    return count

def main():

    # Neighbour joining tree - combined with edit distance:

    n_a = 2
    n_b = 2
    
    # Create random matrix
    pairwise_alignment_distance = create_random_pairwise_alignment_distance_matrix(n_a,n_b)
    
    # Create save copy of matrix
    tmp = copy(pairwise_alignment_distance)        
    
    events = edit_distance_nj(pairwise_alignment_distance, n_a, n_b)
    
    print(tmp)
    print(events)


########################################### MAIN #########################################   
        
if __name__=="__main__":
    
    main() 