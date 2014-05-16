import itertools, logging, os, random, shutil, subprocess, sys, tempfile
import scipy

from repeat.paths import *
from . import align,ensembl_species_info

sys.path.append(os.path.join(ROOT,'lib/python3.2.2/'))
from ete2 import Tree

LOGLEVEL_VERBOSE = 8
logger = logging.getLogger(__name__)
logging.addLevelName(LOGLEVEL_VERBOSE, "VERBOSE") # for parser debug output


################################# TREE TRANSVERSAL #######################################

SPECIES_TO_ID = {'Saccaromyces_cerevisiae': 'Y', 'Oreochromis_niloticus': 'ENSONIP', 'Bos_taurus': 'ENSBTAP', 'Gorilla_gorilla': 'ENSGGOP0', 'Pan_troglodytes': 'ENSPTRP0', 'Homo_sapiens': 'ENSP0', 'Mus_musculus': 'ENSMUSP0', 'Loxodonta_africana': 'ENSLAFP000'}
# Also: 'Saccaromyces_cerevisiae': 'C'  

def get_paralogs(tree_newick, species = 'Loxodonta_africana'):

    '''
    Get all paralogs in <tree_newick> corresponding to <species>.
  
    Parameters:
    <tree_newick> can either be a newick tree string, or the path to a file in newick format.
    
    Return:
    List of identifiers
    '''

    if result_species in SPECIES_TO_ID:
        result_identifier = SPECIES_TO_ID[result_species]
    else:
        raise ValueError('result_species: {0} not in SPECIES_TO_ID'.format(result_species))
        return None

    # Load tree
    tree = Tree(tree_newick, format = 1)   

    paralogs = []
    for iNode in tree.iter_leaves():
        if iNode.name.startswith(result_identifier):
            paralogs.append(iNode.name)
    return paralogs
        

def get_ortholog_pairs(tree_newick, lSpecies):

    '''
    Current algorithm: 
    - Find all paralogs for both species in tree_newick
    - Find for all paralogs in species 1 the corresponding ortholog in species 0, if existant.
    - Check for all paralogs in species 0, that do not have an ortholog yet, whether they should be part of the orthologous pair
     (I.e. because a duplication event has occured after the speciation.)
    '''
    
    try:
        lSpecies_identifier = [SPECIES_TO_ID[iSpecies] for iSpecies in lSpecies]
    except:
        raise ValueError('a member of lSpecies: {0} not in SPECIES_TO_ID'.format(lSpecies))    
        return None
    
    # Load tree
    tree = Tree(tree_newick, format = 1)    
    
    # Find all paralogs of both species in lSpecies in tree:
    paralogs_0 = []
    paralogs_1 = []
    for iNode in tree.iter_leaves():
        if iNode.name.startswith(lSpecies_identifier[0]):
            paralogs_0.append(iNode)
        elif iNode.name.startswith(lSpecies_identifier[1]):
            paralogs_1.append(iNode)
    
    if paralogs_0 == [] or paralogs_1 == []:
        return None
        
    unique_orthology = False    
    pairs = {}
    for i,iP_0 in enumerate(paralogs_0):
        # Detect orthologs to iP_0 in paralogs_1 
        orthologs = []
        for iP_1 in paralogs_1:
            if tree.get_common_ancestor(iP_0,iP_1).D == 'N':
                orthologs.append(iP_1)
        if len(orthologs) > 0:
            pairs[iP_0] = orthologs
            for iO in orthologs:
                paralogs_1.remove(iO)
        if len(paralogs_1) == 0:
            paralogs_0_remainder = paralogs_0[i+1:]
            break
    else:
        # All established orthology relationships in pairs_0 are unique.
        # The remaining paralogs in paralogs_1 do not have a ortholog in species_0
        unique_orthology = True
        
    # Transform pairs_0 to list of lists:
    ortholog_pairs = [[[iP_0],orthologs] for iP_0, orthologs in pairs.items()]
        
    if not unique_orthology:    
        # Check for the remaining paralogs <paralogs_0_remainder> whether they are ortholog to any member of pairs_0:
        for iP_0_r in paralogs_0_remainder:
            for lP_0, lP_1 in ortholog_pairs:
                if tree.get_common_ancestor(iP_0_r,lP_1[0]).D == 'N':
                    lP_0.append(iP_0_r)
    
    ortholog_pairs = [[[iP.name for iP in lP_0],[iP.name for iP in lP_1]] for lP_0,lP_1 in ortholog_pairs]
    
    return  ortholog_pairs 
    
def get_closest_ortholog(gene_tree_newick, species_tree_newick, start_ID):

    start_species = ensembl_species_info.translate('ensembl_ID', 'official', start_ID)
    species_tree = Tree(species_tree_newick, format = 1)  
    lO = get_orthologs(gene_tree_newick, start_ID)
    # clade_name,clade_distance
    lClade = [get_clade_info(species_tree_newick,['Homo_sapiens',ensembl_species_info.translate('ensembl_ID', 'official', iO)]) for iO in lO]
    lClade_Index = [i[1] for i in lClade]
    iO_index_closest = lClade_Index.index(min(lClade_Index))
    
    # [closest_ortholog, clade_name, clade_distance]
    return lO[iO_index_closest], lClade[iO_index_closest]

    

def get_orthologs(tree_newick, start_ID, lSpecies = None):
    
    ''' 
    Get all orthologs to <start_ID> on tree <tree_newick>.
    If <lSpecies> is defined, return only orthologs belonging to these species.
    Current algorithm: Search all leaves belonging to <lSpecies>. 
    If the ancestral node is a speciation, save the leaves as orthologs.
    
    
    Parameters:
    <tree_newick> can either be a newick tree string, or the path to a file in newick format.
    <start_ID> is the node.name of the start entry
    <lSpecies>  ['Loxodonta_africana'] is the species of the output entry
    
    Output:
    List of identifiers
    '''
    
    # Load tree
    tree = Tree(tree_newick, format = 1)    
    
    # Get node corresponding to start_id. Assume that start_id is a unique identifier.
    start_node = tree.get_leaves_by_name(start_ID)[0]
    
    orthologs = []
    for iNode in tree.iter_leaves():
        if iNode != start_node and tree.get_common_ancestor(start_node,iNode).D == 'N':
            orthologs.append(iNode.name)
        
    if lSpecies:       
        try:
                lID = [SPECIES_TO_ID[i] for i in lSpecies]       
        except:
            raise ValueError('species: {0} not in SPECIES_TO_ID'.format(lSpecies))
            return None
        
        orthologs = [iO for iO in orthologs if any(iO.startswith(iID) for iID in lID)]          
                   
    return orthologs
    
def get_within_species_paralogs(tree_newick, species):

    ''' Find all within <species> paralogs
        and return as a list        
    '''

    # Load tree
    tree = Tree(tree_newick, format = 1)    
    
    # Return within <species> paralogs      
    return [iLeaf.name for iLeaf in tree.iter_leaves() if ensembl_species_info.translate('ensembl_ID', 'official', iLeaf.name) == species]



def get_inparalogs(tree_newick, species):

    ''' Find all inparalogs, i.e. genes within <species>, that have arisen through duplications
        after the last speciation.
        Return all groups of inparalogs as a list of lists.
        
    '''
    
    # Load tree
    tree = Tree(tree_newick, format = 1)    
    
    # get_within_species_paralogs
    paralogs = [iLeaf for iLeaf in tree.iter_leaves() if ensembl_species_info.translate('ensembl_ID', 'official', iLeaf.name) == species]
    
    lInparalogs = []
    
    while paralogs:
        inparalogs = [paralogs.pop()]
        for iP in paralogs:
            ancestor = tree.get_common_ancestor(inparalogs[0],iP)
            current_parent = iP.up
            while current_parent != ancestor:
                if current_parent.D == 'N':
                    break
                current_parent = current_parent.up
            else:
               inparalogs.append(iP)
               paralogs.remove(iP)
        inparalogs = [iP.name for iP in inparalogs]
        lInparalogs.append(inparalogs)
        
    return lInparalogs

    

def get_orthologs_before_first_duplication(tree_newick, start_id, result_species = 'Loxodonta_africana'):
    
    ''' 
    Get all orthologs to <start_id> in species <result_species> on tree <tree_newick>, that
    have split more recently than the most recent duplication.
    Current algorithm: Search all leaves below the closest duplication node. Add any of them
    belonging to <result_species> to the list of orthologs.
    
    Parameters:
    <tree_newick> can either be a newick tree string, or the path to a file in newick format.
    <start_id> is the node.name of the start entry
    <result_species> is the species of the output entry
    
    Output:
    List of identifiers
    '''
    
    if result_species in SPECIES_TO_ID:
        result_identifier = SPECIES_TO_ID[result_species]
    else:
        raise Error('result_species: {0} not in SPECIES_TO_ID'.format(result_species))
        return None
    
    # Load tree
    tree = Tree(tree_newick, format = 1)    
    
    # Get node corresponding to start_id. Assume that start_id is a unique identifier.
    start_node = tree.get_leaves_by_name(start_id)[0]
    
    # Search first duplication node above <start_node>
    clostest_duplication_node = start_node
    while clostest_duplication_node.D != 'Y':
        clostest_duplication_node = clostest_duplication_node.up
     
    orthologs = []
    for iNode in clostest_duplication_node.iter_leaves():
        if iNode.name.startswith(result_identifier):
            orthologs.append(iNode.name)
    return orthologs

###################### Functions for ORDER CONSERVATION ##################################


def consecutive_ranges(data):
    ''' find range of consecutive integers in <data>. E.g.
        data = [0,1,4,5,6,10] returns
        [[0,1],[4,6],[10,10]]
        CURRENTLY NOT USED
        '''
    
    from operator import itemgetter    
    data.sort()  
    # From http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    result = []
    for k, g in itertools.groupby(enumerate(data), lambda x:x[0]-x[1]):
        range = list(map(itemgetter(1), g))
        result.append([range[0],range[-1]])
    return result
    

def get_clade_info(tree_newick,lNode_names):

    ''' return distance to and name of MRCA (clade) node'''
    
    if lNode_names[0] == lNode_names[1]:
        return lNode_names[0],-1
    
    tree = Tree(tree_newick, format = 1)  
    try:
        ancestor = tree.get_common_ancestor(lNode_names)
    except:
        raise ValueError('Perhaps lNode_names: {0} not in tree:\n{1}'.format(lNode_names,str(tree)))
        return None 
    
    distance = 0
    current_node = tree.get_leaves_by_name(name=lNode_names[0])[0]
    while current_node.up != ancestor:
        distance += 1
        current_node = current_node.up

    return ancestor.name,distance

def get_node_order(tree_newick,start_node_name):

    ''' return list of nodes in <tree_newick> ordered by their distance to <main_node> '''
    
    tree = Tree(tree_newick, format = 1)
    #tree.unroot()
    
    try:
        start_node = tree.get_leaves_by_name(name=start_node_name)[0]
    except:
        raise ValueError('start_node: {0} not in tree:\n{1}'.format(start_node,str(tree)))
        return None   
        
    order = []    
    current_node = start_node
    while not current_node.is_root():
        order.append({current_node.up.name : list(current_node.get_sisters()[0].get_leaf_names())})
        current_node = current_node.up
        
    return order
    



################################ TREE RECONSTRUCTION #####################################

def reconstruct_tree(homologs, method = 'PhyML', sequence_type = 'aa', result_file = None):

    '''
    Multiple sequence align <sequences>. If <tandem_repeat_annotations> are given, use <prograph>
    and provide them as added information in t-reks format.
    
    Parameters: Dict of sequences and identifiers
    e.g. {'ENSP00012': 'ABCABAC', 'ENSP00013': 'BCABCA'}   '''
    
    if len(homologs) == 0:
        raise ValueError('Dict <homologs> is empty.')
        return None
    
    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    logger.debug("Reconstruct_tree: Created temp directory: %s", working_dir)

    # Save MSA to file in temp directory
    alignment_file = os.path.join(working_dir, 'input_alignment.phylip')
    align.save_msa_phylip(homologs, alignment_file)
     
    if method == 'PhyML':
        # Run PhyML
        # See http://www.atgc-montpellier.fr/download/papers/phyml_manual_2009.pdf for choice of options.
        # The PhyML result is in stdout?
        phyML_path = os.path.join(EXECROOT, 'PhyML_3.0_linux64')
        tree_file = os.path.join(working_dir,'input_alignment.phylip_phyml_tree.txt')
        p = subprocess.Popen([phyML_path, "-i", alignment_file, "-d", sequence_type, "-q", "-s NNI", "2>/dev/null"])
        p.wait()
        
        if result_file:
            shutil.copyfile(tree_file, result_file)
        
        #from Bio import Phylo
        #tree = Phylo.read(tree_file, "newick")
        
        with open(tree_file,'r') as th:
            tree = th.readline().rstrip()
        logger.debug(tree)
        
        return tree

    
######################################### MAIN ###########################################

def main():
    start_id = 'ENSPCAP00000006957'
    tree_newick = '/cluster/home/infk/anmaria/ex.nhx'
    get_orthologs(tree_newick, start_id, result_species = 'Loxodonta_africana')