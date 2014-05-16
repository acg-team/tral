# Copyright (C) 2012-2013 by Elke Schaper (elke@inf.ethz.ch) 

"""
    Analyse tandem repeat phylogenies
    of homolog tandem repeats in genes within a gene tree
    
""" 

import itertools,logging,os,pickle,re
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import scipy as sc
from scipy.stats.distributions import binom as binom
from scipy.stats import hypergeom

import logging
logger = logging.getLogger(__name__)
#logging.basicConfig(level=logging.INFO)

from . import align, ensembl_IO, gene_tree_info, tree_io, tree_operations, ensembl_species_info, repeat_tree_info
from repeat.paths import *
from ..hmm import hmm, sequence, viterbi
from ..repeat import repeat_info, repeat_realign, repeat_io, repeat_exon
    
dTranslate = {  'ortholog_present': "[{'Ortholog': None}]",
           'furthest_TR_n_4': "[{'Ortholog_n': ['min', 4]}]", 
           'furthest_TR_p': "[{'dn': ['max', 0], 'Kendall': ['min', 1], 'd(max_n,nBC)': ['max', 0]}]",
           'furthest_TR_OG_d1_K1' : "[{'Kendall': ['min', 1], 'n1': ['min', 4], 'd(max_n,nBC)': ['max', 1], 'nBC': ['min', 4]}]",
           'closest_TR_none_perfect': "[{'dn': ['min', 1]}, {'Kendall': ['less', 1]}, {'d(max_n,nBC)': ['min', 1]}]", 
           'closest_TR_completely_diverged': "[{'n1': ['min', 4], 'min_split': ['max', 1]}]",
           'closest_TR_strongly_diverged_1': "[{'n1': ['min', 4], 'min_split': ['max', 2]}]",
           'closest_TR_strongly_diverged_2': "[{'n1': ['min', 4], 'min_split': ['max', 2]}, {'dn': ['min', 4], 'n1': ['min', 4]}]",
           'closest_max_OnD_1': "[{'Ortholog_nD': ['max', 1]}]",
           'closest_max_OnD_0_8': "[{'Ortholog_nD': ['max', 0.8]}]",
           'closest_min_dn_4': "[{'dn': ['min', 4]}]",
           'closest_difference_units_4': "[{'Ortholog_dn': ['min', 4]}]",
           'closest_difference_units_5': "[{'Ortholog_dn': ['min', 5]}]",
           'closest_TR_partially_diverged': "[{'dn': ['min', 4], 'd(max_n,parsimony)': ['min', 4]}]",
           }

dPFAM_Name = {'PF00400': 'WD40 repeat',
        'PF00023': 'Ankyrin repeat',
        'PF07646': 'Kelch motif 2',
        'PF01344': 'Kelch motif 1',
        'PF00515': 'Tetratricopeptide 1',
        'PF00514': 'Armadillo repeat',
        'PF00090': 'Thrombospondin',
        'PF07719': 'Tetratricopeptide 2',
        'PF07679': 'Immunoglobulin I-set',
        'PF07686': 'Immunoglobulin V-set',
        'PF00076': 'RNA recognition motif',
        'PF00045': 'Hemopexin',
        'PF00520': 'TM ion channel',
        'PF00415': 'RCC1',
        'PF02210': 'Laminin',
        'PF00595': 'PDZ domain',
        'PF00013': 'KH domain',
        'PF00041': 'Fibronectin type III',
        'PF00096': 'Zinc finger',
        'PF12171': 'JAZ Zinc finger',
        'PF00084': 'Sushi domain/SRC repeat',
        'PF10473': 'Cenp-F/LEK1 LRR',
        'PF01391': 'Collagen',
        't-reks' : 'T-REKS',
        'hhrepid': 'HHREPID',
        'xstream': 'XSTREAM',
        'PF00628': 'PHD-finger',
        'PF00059': 'Lectin C-type domain',
        'PF00412': 'LIM domain',
        'PF07653': 'SH3 domain 2',
        'PF02493': 'MORN repeat',
        'PF02985': 'HEAT repeat',
        'PF00028': 'Cadherin',
        'PF06758': '\emph{unknown function}',
        'PF00024': 'PAN domain',
        'PF03516': 'Filaggrin',
        'PF07142': '\emph{unknown function}',
        'PF07974': 'EGF-like domain 2',
        'PF07926': 'TPR/MLP1/MLP2-like protein',
        'PF00418': 'tubulin-binding repeat',
        'PF00098': 'Zinc knuckle',
        'PF00307': 'Calponin homology (CH) domain',
        'PF00432': 'Prenyltransferase and squalene oxidase repeat',
        'PF00637': 'Clathrin repeat',
        'PF00036': 'EF hand',
        'PF00560': 'Leucine rich repeat',
        'PF07645': 'Calcium-binding EGF domain',
        'PF01851': 'Proteasome/cyclosome repeat',
        'PF00435': 'Spectrin repeat',
        'PF00191': 'Annexin',
        'PF00020': 'TNFR/NGFR cysteine-rich region',
        'PF05386': 'TEP1 N-terminal domain',
        'PF00240': 'Ubiquitin',
        'PF11976': 'Ubiquitin-2 like Rad60 SUMO-like',
        'PF08065': 'K167R (NUC007) repeat',
        'PF00612': 'IQ calmodulin-binding motif',
        'PF00018': 'SH3 domain',
        'PF08205': 'CD80-like C2-set immunoglobulin domain',
        'PF02469': 'Fasciclin domain',
        'PF07004': 'Sperm-tail PG-rich repeat',
        'PF00093': 'Von Willebrand factor type C domain',
        'PF00047': 'Immunoglobulin domain',
        'PF01826': 'Trypsin Inhibitor like cysteine rich domain',
        'PF07648': 'Kazal-type serine protease inhibitor domain',
        'PF08913': 'Vinculin Binding Site',
        'PF02185': 'Hr1 repeat',
        'PF10629': '\emph{unknown function}',
        'PF04680': 'Opioid growth factor receptor',
        'PF01390': 'SEA domain',
        'PF02818': 'PPAK motif',
        'PF07720': 'Tetratricopeptide repeat 3',
        'PF11627': 'Nuclear factor hnRNPA1',
        'PF01437': 'Plexin repeat',
        'PF00008': 'EGF-like domain',
        'PF00053': 'Laminin EGF-like',
        'PF00094': 'von Willebrand factor type D domain',
        }


def create_tree(tandem_repeat_ID, my_TR_data, my_gene_tree, lID, result_dir, key = None, create_draw_tree_files = None):
    
    ''' Find all viterbi paths for <lID> in <my_TR_data>.
        Translate to sequence data with sequence in <my_gene_tree>.
        
        Realign MSAs with Mafft.
        Build phylogeny with PhyML
        Create .pdf of phylogeny with scriptree.
        
        Save resulting files to <result_dir>, with <tandem_repeat_ID> as first identifier,
        and an ID in <lID> as second identifier.
        
        Return the aligned MSA, and the phylogeny in newick format.
        
         '''
    
    if not key:
        key = "{0}_{1}".format(str(tandem_repeat_ID),lID[0])
    
    # Find all Viterbi paths
    lViterbi_path = [my_TR_data['homologs'][ensembl_ID].viterbi_path for ensembl_ID in lID if type(my_TR_data['homologs'][ensembl_ID]) != dict]
    if len(lViterbi_path) <= 1:
        return None, None
    
    # Reconstruct the repeat unit MSA from the Viterbi Paths of homologous tandem repat prediction.
    lSequence = [str(my_gene_tree['leafs'][ensembl_ID]['sequences']['aa'].seq).replace("*","").replace("X","") for ensembl_ID in lID]
    
    lMSA = viterbi.hmm_path_to_maximal_complete_tandem_repeat_units(lSequence, lViterbi_path, my_TR_data['lD'], alpha = 0.6)
    
    # Realign the repeat unit MSA using MAFFT
    repeat_units = {}
    annotations = {}
    for ensembl_ID,iMSA in zip(lID,lMSA):
        for i,iUnit in enumerate(iMSA):
            repeat_units['{0}_{1}'.format(ensembl_ID, str(i))] = iUnit
            species = ensembl_species_info.translate('ensembl_ID','official', ensembl_ID)
            annotations['{0}_{1}'.format(ensembl_ID, str(i))] = 'Species {{{0}}} ensembl_ID {{{1}}} RepeatUnit {{{2}}}'.format(species, ensembl_ID, str(i))
    if len(repeat_units) <= 1:
        return None, None
    
    aligned_homologs = align.align(sequences = repeat_units, aligner = 'MAFFT')
    if len(aligned_homologs) <= 1:
        return None, None
    
    aligned_msas = {ensembl_ID: [aligned_homologs['{0}_{1}'.format(ensembl_ID, str(i))]  for i in range(len(iMSA)) ] for ensembl_ID,iMSA in zip(lID,lMSA)}
    alignment_file = os.path.join(result_dir, "{0}.faa".format(key))
    #repeat_io.save_repeat_fasta(aligned_msas, alignment_file)
    
    # Reconstruct the repeat unit tree using PhyML
    #tree_file = os.path.join(result_dir, "{0}.nhx".format(key))
    tree_file = None
    newick_tree = tree_operations.reconstruct_tree(aligned_homologs, result_file=tree_file)
    
    if create_draw_tree_files:
        # Create scripttree files to later create script trees.
        tree_draw_file = os.path.join(result_dir, "{0}".format(key))
        #tree_io.create_scriptree_files(result_file_stump = os.path.join(result_dir, "{}_{}_{}".format(my_TR_ID,ensembl_ID,iO)), tree = my_TR_data['homologous_pairs'][ensembl_ID][iO], annotations = None)
        tree_io.create_scriptree_files(tree_draw_file, tree_nhx = newick_tree, annotations = annotations)
    
    return aligned_msas, newick_tree
    
 
def get_pairwise_ortholog_repeat_unit_tree(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, lSpecies):
    
    '''
        Find all orthologous pairs of lSpecies in the gene_tree corresponding to my_TR.
        For each pair:
            - Reconstruct the repeat unit MSA from the Viterbi Paths of homologous tandem repat prediction.
            - Realign the repeat unit MSA using MAFFT
            - Reconstruct the repeat unit tree using PhyML
            - Save reconstructed tree as file, and as image (scripttree)
            
        Parameters:
        <lSpecies> is a list of two species
    '''
    
    ortholog_dir = os.path.join(tandem_repeat_dir, "orthologous_pairs", "_".join(lSpecies))
    if not os.path.exists(ortholog_dir):
        os.makedirs(ortholog_dir)
    
    my_TR_data, my_gene_tree = gene_tree_info.load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)
    results = {'tandem_repeat': tandem_repeat_ID, 'gene_tree': my_TR_data['gene_tree'], 'species': lSpecies}
    
    # Find all orthologous pairs in the gene_tree corresponding to my_TR.
    orthologous_pairs = tree_operations.get_ortholog_pairs(my_gene_tree['tree_nhx'], lSpecies)
    if not orthologous_pairs:
        return None
    results['orthologous_pairs'] = {}
    
    for iO in orthologous_pairs:
        # Merge IDs
        lO = list(itertools.chain(*iO))
        results['orthologous_pairs'][lO[0]] = {'all_identifier': lO}
        
        ### missing
        aligned_msas, newick_tree = create_tree(tandem_repeat_ID, my_TR_data, my_gene_tree, lO, ortholog_dir)
        results['orthologous_pairs'][lO[0]]['alignment'] = aligned_msas
        results['orthologous_pairs'][lO[0]]['tree'] = newick_tree

            
    with open(os.path.join(ortholog_dir, str(tandem_repeat_ID)+".pickle"), 'wb') as result_handle:
        pickle.dump(results, result_handle)
        

def order_conservation_pairwise_ortholog_repeat_unit_tree(tandem_repeat_ID, tandem_repeat_dir, lSpecies, result_file):
    
    ortholog_dir = os.path.join(tandem_repeat_dir, "orthologous_pairs", "_".join(lSpecies))
    
    try:
        with open(os.path.join(ortholog_dir,str(tandem_repeat_ID) + '.pickle'),'rb') as rh:
            orthologous_data = pickle.load(rh)
        orthologous_pairs = orthologous_data['orthologous_pairs']
    except:
        return None
    
    my_TR_data = load_pickles(tandem_repeat_ID, tandem_repeat_dir)
    
    nPermutations = 1000
    with open(result_file, 'a') as fh:
        for ensembl_ID,iO in orthologous_pairs.items():
            if iO['tree'] == '':
                continue
            # Use only two identifier for each species:
            lSpeciesID = [SPECIES_TO_ID[iSpecies] for iSpecies in lSpecies]
            for iID,jID in itertools.combinations(iO['all_identifier'],2):
                lID = [iID,jID]
                if not ((lID[0].startswith(lSpeciesID[0]) and lID[1].startswith(lSpeciesID[1])) or (lID[0].startswith(lSpeciesID[1]) and lID[1].startswith(lSpeciesID[0]))):
                    continue
                n0 = len(orthologous_pairs[ensembl_ID]['alignment'][lID[0]])   
                n1 = len(orthologous_pairs[ensembl_ID]['alignment'][lID[1]])
                if n0 == 0 or n1 == 0:
                    continue           
                leaves = {id+"_"+str(i):[id,i] for id in lID for i in range(max(n0,n1))}    
                try:
                    result, tree = tree_operations.order_conservation(tree_newick = iO['tree'], leaves=leaves, nPermutation = nPermutations)    
                except:
                    print("leaves = {0}".format(leaves))
                    print("tree_newick = '{0}'".format(iO['tree']))
                    print("lID = {0}".format(lID))
                    print("tandem_repeat_ID = {0}".format(tandem_repeat_ID))
                # ['tandem_repeat_ID','ortholog_group_ID', 'ensembl_ID_0', 'ensembl_ID_1', 'lD', 'n_0', 'n_1', 'nSpeciation', 'Kendall', 'Spearman', 'nPermutation', 'pValue_Kendall', 'pValue_Spearman']
                #
                fh.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}\n".format(tandem_repeat_ID, 
                        ensembl_ID, lID[0], lID[1], my_TR_data['lD'], n0, n1, 
                        result['nSpeciation'], result['kendall_biocluster'], result['spearman_biocluster'],
                        nPermutations, result['reshuffle']['kendall_biocluster'], result['reshuffle']['spearman_biocluster']))


def phylogeny_all_orthologs(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, result_dir, species_main = 'Homo_sapiens', lSpecies = ['Gorilla_gorilla', 'Mus_musculus']):
    
    '''
    Get orthologs and within species paralogs for <species_main>. 
    If <lSpecies> is defined, search only for orthologs to any species in <lSpecies>.
    
    Assemble all PhyML tandem repeat unit trees for all pairwise orthologs, and pairs of within species paralogs.
    These trees as repeat_tree_info.Repeat_tree. 
    '''

    result_file = os.path.join(result_dir,str(tandem_repeat_ID)+'.pickle')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    my_TR_data, my_gene_tree = gene_tree_info.load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)
    
    if not 'homologs' in my_TR_data:
        logging.error('No homologs were found to this TR.')
        return None
    
    my_TR_data['homologous_pairs'] = {}
    
    lTR_containing = [ensembl_ID for ensembl_ID, repeat_homolog in my_TR_data['homologs'].items() if not type(repeat_homolog) == dict and repeat_homolog.nD >= 0]
    lInparalogs = tree_operations.get_inparalogs(my_gene_tree['tree_nhx'], species_main)
    lInparalogs = filter(lambda x: x != [], [[i for i in iInparalogs if i in lTR_containing] for iInparalogs in lInparalogs])
    lWithin_Species_Paralogs = tree_operations.get_within_species_paralogs(my_gene_tree['tree_nhx'], species_main)
    #lWithin_Species_Out_Paralogs = [i for i in lWithin_Species_Paralogs if i not in lInparalogs]
        
    species_tree = os.path.join(DATAROOT2, 'Compara','species_tree_ncbi.nhx')
    #species_order = tree_operations.get_node_order(species_tree,species_main)
    
    for iInCount, iInparalogs in enumerate(lInparalogs):
        orthologs = tree_operations.get_orthologs(my_gene_tree['tree_nhx'],iInparalogs[0],lSpecies)
        orthologs = [i for i in orthologs if i in lTR_containing]
        for iP, iO in itertools.product(iInparalogs, orthologs):
            ortholog_species = ensembl_species_info.translate('ensembl_ID', 'official', iO)
            clade_name,clade_distance = tree_operations.get_clade_info(species_tree,lNode_names = [species_main,ortholog_species])
            lAligned_msa, newick_tree = create_tree(tandem_repeat_ID, my_TR_data, my_gene_tree, lID = [iP,iO], result_dir = result_dir, key = '{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(tandem_repeat_ID,iInCount,clade_distance,clade_name,ortholog_species,iP, iO))
            if not newick_tree:
                continue
            if not iP in my_TR_data['homologous_pairs']:
                my_TR_data['homologous_pairs'][iP] = {}
            iTree = create_repeat_unit_tree(newick_tree,clade_name,clade_distance,lAligned_msa)
            if not iTree == None:
                my_TR_data['homologous_pairs'][iP][iO] = iTree
    for iP, iO in itertools.product(lWithin_Species_Paralogs, lWithin_Species_Paralogs):
        if iP == iO:
            continue
        ortholog_species = species_main
        clade_name,clade_distance = species_main,0
        lAligned_msa, newick_tree = create_tree(tandem_repeat_ID, my_TR_data, my_gene_tree, lID = [iP,iO], result_dir = result_dir, key = '{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(tandem_repeat_ID,iInCount,clade_distance,clade_name,ortholog_species,iP, iO))
        if not newick_tree:
            continue
        if not iP in my_TR_data['homologous_pairs']:
            my_TR_data['homologous_pairs'][iP] = {}
        iTree = create_repeat_unit_tree(newick_tree,clade_name,clade_distance,lAligned_msa)
        if not iTree == None:
            my_TR_data['homologous_pairs'][iP][iO] = iTree
    
    with open(result_file, 'wb') as result_handle:
        pickle.dump(my_TR_data, result_handle)
  

def create_repeat_unit_tree(newick_tree,clade_name,clade_distance,msa):
    
    ''' Create a Repeat_Tree instance. 
        Calculate cherries and bisample cherries, and a premature order conservation test
        At everything after <newick_tree> as tree attributes. LATER: Do not hardcode these attributes.
    '''
    
    if newick_tree == '':
        return None
    
    iTree = repeat_tree_info.Repeat_Tree(newick_tree)    
    iTree.unroot()
    iTree.set_Ensembl_ID_and_index()
    if len(iTree.lID) < 2:
        return None
    iTree.get_cherries()
    iTree.get_bisample_cherries()
    iTree.get_min_split()
    iTree.order_conservation_test()
    iTree.clade_name = clade_name 
    iTree.clade_distance = clade_distance
    iTree.msa = msa
    
    return iTree

def create_pairwise_gene_tree_data(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, result_file, species_main = 'Homo_sapiens'):
    
    ''' Do all the calculations of pairwise gene tree data, that has otherwise been missing. Can surely be shortened by a lot.''' 
    
    my_TR_data, my_gene_tree = gene_tree_info.load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)
        
    my_TR_data['conservation'] = {}
    
    # Get all species in gene_tree
    lSpecies = [ensembl_species_info.translate("ensembl_ID", "official", ensembl_ID) for ensembl_ID in my_gene_tree['leafs'].keys()]
                
    # Get distances on gene tree
    species_tree = os.path.join(DATAROOT2, 'Compara','species_tree_ncbi.nhx')
    lClade = [tree_operations.get_clade_info(species_tree,lNode_names = [species_main,iSpecies]) for iSpecies in lSpecies] 
    
    my_TR_data['lHMM'] = my_TR_data['lD']
    lHMM = my_TR_data['lHMM']
     
    my_TR_data['conservation']['tr_present_oldest'] = {}
    
    my_TR_data['orthology'] = {}
    
    my_TR_data['conservation']['paralog_wise'] = {}
    lWithin_Species_Paralogs = tree_operations.get_within_species_paralogs(my_gene_tree['tree_nhx'], species_main)
    
    # Calculate mean TR unit length
    lMean = [my_TR_data['homologs'][iP].lD for iP in lWithin_Species_Paralogs if type(my_TR_data['homologs'][iP]) != dict]
    my_TR_data['lMean'] = sum(lMean)/len(lMean)
    
    for iP in lWithin_Species_Paralogs:
        my_TR_data['conservation']['paralog_wise'][iP] = {}
        iP_data = my_TR_data['conservation']['paralog_wise'][iP]
        
        # Calculate divergence
        if type(my_TR_data['homologs'][iP]) != dict and not hasattr(my_TR_data['homologs'][iP], 'divergence'):
            my_TR_data['homologs'][iP].calculate_scores(scoreslist=['phylo_gap001'])
                
        # Ortholog present
        lOrtholog = tree_operations.get_orthologs(my_gene_tree['tree_nhx'],iP)
        lOrtholog_species = [ensembl_species_info.translate("ensembl_ID", "official", ensembl_ID) for ensembl_ID in lOrtholog]
        lClade = [tree_operations.get_clade_info(species_tree,lNode_names = [species_main,iSpecies]) for iSpecies in lOrtholog_species]
        my_TR_data['orthology'][iP] = {'lOrtholog': lOrtholog, 'lClade':lClade}
        
        # Do missing phylogeny calculations.
        for iO in lOrtholog:
            try:
                iTree = my_TR_data['homologous_pairs'][iP][iO]
            except:
                logging.debug("Key combination {}, {} does not exist".format(iP,iO))
                continue
            if iTree == None:
                continue    
            iTree.get_bisample_cherries()
            iTree.get_min_split()
        
        # Exon structure
        iP_data['exon'] = {}
        exon_structure = my_gene_tree['leafs'][iP]['exon']
        my_TR = my_TR_data['homologs'][iP] 
        try:
            nExon, nTR_in_Exon_Max = repeat_exon.get_exon_measures(my_TR, exon_structure)  
        except:
            # my_TR is probably not a tandem repeat object.
            nExon = 0
            nTR_in_Exon_Max = 0
        iP_data['exon'] = {'nExon': nExon, 'nTR_in_Exon_Max': nTR_in_Exon_Max}
        
    
    with open(result_file, 'wb') as result_handle:
        pickle.dump(my_TR_data, result_handle) 


def get_attributes(my_TR_data, ensembl_ID, attributes, ortholog = None):
    
    if ortholog:
        if ortholog in my_TR_data['homologous_pairs'][ensembl_ID]:
            pairwise_tree = my_TR_data['homologous_pairs'][ensembl_ID][ortholog]
        ortholog_TR = my_TR_data['homologs'][ortholog]
                    
    iAttributes = []
    for iA in attributes:
        ######  Human TR ######
        if iA == 'lHMM': 
            iAttributes.append(my_TR_data['lHMM'])
        if iA == 'lMean':
            iAttributes.append(my_TR_data['lMean'])  
        elif iA == 'n':
            try:
                iAttributes.append(my_TR_data['homologs'][ensembl_ID].n)
            except:
                iAttributes.append(1)  
        elif iA == 'nD':
            try:
                iAttributes.append(my_TR_data['homologs'][ensembl_ID].nD)
            except:
                iAttributes.append(1)
        elif iA == 'detection_ID':
            if 'PFAM_ID' in my_TR_data['HMM']:
                iAttributes.append(my_TR_data['HMM']['PFAM_ID'])
            else:
                iAttributes.append(my_TR_data['HMM']['tandem_repeat_seed'].TRD)
        elif iA == 'nExon':
            iAttributes.append(my_TR_data['conservation']['paralog_wise'][ensembl_ID]['exon']['nExon'])
        elif iA == 'nTR_in_Exon_Max':
            iAttributes.append(my_TR_data['conservation']['paralog_wise'][ensembl_ID]['exon']['nTR_in_Exon_Max'])
        elif iA == 'Ensembl_Protein_ID':
            iAttributes.append(ensembl_ID)
        elif iA == 'divergence':
            if not hasattr(my_TR_data['homologs'][ensembl_ID], 'divergence'):
                my_TR_data['homologs'][ensembl_ID].calculate_scores(scoreslist=['phylo_gap001'])
            iAttributes.append(my_TR_data['homologs'][ensembl_ID].divergence['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'])
        elif iA == 'GO':
            # not finished
            infile = open(os.path.join(DATAROOT2,"homo_sapiens_go_protein.tsv"),'r')
            go_terms = ensembl_IO.get_data(infile, dict_ID=('Ensembl Protein ID',ensembl_ID), all=True)
            try:
                iAttributes += [[i['GO Term Accession']] for i in go_terms if 'GO Term Accession' in i]
            except:
                print(go_terms)
        ######  Ortholog TR ###### 
        elif iA == 'Ortholog_lD':
            if type(ortholog_TR) == dict:
                if 'repeat_unit' in ortholog_TR:
                    iAttributes.append(len(ortholog_TR['repeat_unit']))
                else:
                    iAttributes.append(0)
            else:
                iAttributes.append(ortholog_TR.lD)
        elif iA == 'Ortholog_n':
            if type(ortholog_TR) == dict:
                if 'repeat_unit' in ortholog_TR:
                    iAttributes.append(1)
                else:
                    iAttributes.append(0)
            else:
                iAttributes.append(ortholog_TR.n)
        elif iA == 'Ortholog_nD':
            if type(ortholog_TR) == dict:
                if 'repeat_unit' in ortholog_TR:
                    iAttributes.append(1)
                else:
                    iAttributes.append(0)
            else:
                iAttributes.append(ortholog_TR.nD)
        elif iA == 'Ortholog_dn':
            if type(ortholog_TR) == dict:
                if 'repeat_unit' in ortholog_TR:
                    iAttributes.append(abs(1-iAttributes.append(my_TR_data['homologs'][ensembl_ID].n)))
                else:
                    iAttributes.append(iAttributes.append(my_TR_data['homologs'][ensembl_ID].n))
            else:
                iAttributes.append(abs(ortholog_TR.n-iAttributes.append(my_TR_data['homologs'][ensembl_ID].n)))
        elif iA == 'Ortholog_Ensembl_Protein_ID':
            iAttributes.append(ortholog)               
        ###### Pairwise Measures ###### ['n0', 'n1', 'dn', 'min_split', 'nC','nCB','p_nCB','dnBC_nC','Kendall', 'p_Kendall']
        elif iA == 'n0':
            iAttributes.append(pairwise_tree.n[ensembl_ID])
        elif iA == 'n1':
            iAttributes.append(pairwise_tree.n[ortholog])
        elif iA == 'dn':
            iAttributes.append(abs(pairwise_tree.n[ensembl_ID]-pairwise_tree.n[ortholog]))
        elif iA == 'min_split':
            iAttributes.append(pairwise_tree.min_split)
        elif iA == 'nC':
            iAttributes.append(pairwise_tree.nC)
        elif iA == 'nBC':
            iAttributes.append(pairwise_tree.nBC)
        elif iA == 'p_nBC':
            iAttributes.append(pairwise_tree.pValue['nBC'])
        elif iA == 'dnBC_nC':
            iAttributes.append(pairwise_tree.nC - pairwise_tree.nBC)
        elif iA == 'Kendall':
            if 'Kendall' in pairwise_tree.order_conservation and pairwise_tree.order_conservation['Kendall']:
                iAttributes.append(pairwise_tree.order_conservation['Kendall']['test_statistic'])
            else:
                iAttributes.append('NA')
        elif iA == 'p_Kendall':
            if 'Kendall' in pairwise_tree.order_conservation and pairwise_tree.order_conservation['Kendall'] and 'pValue' in pairwise_tree.order_conservation['Kendall']:
                iAttributes.append(pairwise_tree.order_conservation['Kendall']['pValue'])
            else:
                iAttributes.append('NA')
    return iAttributes


def merge_TR_homology_data(tandem_repeat_dir, result_file_name = 'all/all.pickle'):
    
    ''' Collect basic TR data from all singe TR files and save as pickle.'''
    
    result_file = os.path.join(tandem_repeat_dir,result_file_name)
    
    attributes_Human = ['lHMM', 'lMean', 'n', 'nD', 'detection_ID','nExon', 'nTR_in_Exon_Max','Ensembl_Protein_ID','divergence']
    attributes_Ortholog = ['Ortholog_lD', 'Ortholog_n', 'Ortholog_nD', 'Ortholog_Ensembl_Protein_ID','Ortholog_dn']
    attributes_pairwise = ['n0', 'n1', 'dn', 'min_split', 'nC','nBC','p_nBC','dnBC_nC','Kendall', 'p_Kendall']
    attributes_All = attributes_Human + attributes_Ortholog + attributes_pairwise
    all_TR = {'all_data': {}}
    
    for my_TR_ID, my_TR_data in iterate_tandem_repeat_data(tandem_repeat_dir):
        all_TR['all_data'][my_TR_ID] = {}
        for ensembl_ID in my_TR_data['homologous_pairs'].keys():
            
            all_TR['all_data'][my_TR_ID][ensembl_ID] = {}
            
            for iO in my_TR_data['orthology'][ensembl_ID]['lOrtholog']:
                #print("{} {}".format(ensembl_ID,iO))
                if iO not in my_TR_data['homologs']:
                    print(iO)
                    continue
                if iO in my_TR_data['homologous_pairs'][ensembl_ID]:
                    # Save the following features:
                    # n0, n1, nC, nBC, Kendall, gene_present_oldest,
                    iAttributes =  get_attributes(my_TR_data, ensembl_ID, attributes_All, ortholog = iO)
                else:
                    # No tree was build
                    iAttributes = get_attributes(my_TR_data, ensembl_ID, attributes_Human + attributes_Ortholog, ortholog = iO)
                    iAttributes = iAttributes + ["NA"]*len(attributes_pairwise)
                all_TR['all_data'][my_TR_ID][ensembl_ID][iO] = {name:attr for name,attr in zip(attributes_All, iAttributes)}
    
    with open(result_file, 'wb') as rfh:
        pickle.dump(all_TR,rfh)
    
    return all_TR


def oldest_conserved(data_file, result_file, species_main = 'Homo_sapiens'):
    
    with open(data_file, 'rb') as dfh:
        all_TR = pickle.load(dfh)
             
    all_TR['merged_data'] = {}
    nTotal = len(all_TR['all_data'])
    for i,tmp in enumerate(all_TR['all_data'].items()):
        TR_ID, data_TR = tmp
        print("{}: {}/{}".format(TR_ID, i, nTotal))
        
        all_TR['merged_data'][TR_ID] = {}
        for iP, data_iP in data_TR.items():
            all_TR['merged_data'][TR_ID][iP] = {}
            iParalog_Data = all_TR['all_data'][TR_ID][iP]
            iMerged_Data = all_TR['merged_data'][TR_ID][iP]
            
            ### Most distant observation
            attributes_conservation = []
            # Furthest clade ortholog present
            attributes_conservation.append([{'Ortholog': None}])
            # Furthest clade (Ortholog_nD > 0.8, Ortholog_nD>3.2)
            attributes_conservation.append([{'Ortholog_nD': ['min', 0.8]}])
            attributes_conservation.append([{'Ortholog_nD': ['min', 3.2]}])
            attributes_conservation.append([{'Ortholog_n': ['min', 4]}])
            # Furthest clade perfectly conserved
            attributes_conservation.append([{'dn': ['max', 0], 'Kendall': ['min', 1], 'd(max_n,nBC)': ['max', 0]}])
            # Furthest clade strongly conserved
            attributes_conservation.append([{'n1': ['min', 4], 'd(max_n,nBC)': ['max', 1], 'Kendall': ['min', 1], 'nBC': ['min', 4]}])
            for i in attributes_conservation:
                clade_distance = most_extreme_TR(paralog_Data = iParalog_Data, species_main = species_main, lConstraint = i)
                iMerged_Data[str(i)] = clade_distance
            
            ### Closest observation
            attributes_divergence = []
            # Closest clade completely diverged
            attributes_divergence.append([{'n1': ['min', 4], 'min_split': ['max', 1]}])
            # Closest clade strongly diverged (1)
            attributes_divergence.append([{'n1': ['min', 4], 'min_split': ['max', 2]}])
            # Closest clade strongly diverged (2)
            attributes_divergence.append([{'n1': ['min', 4], 'min_split': ['max', 2]}, {'n1': ['min', 4], 'dn': ['min', 4]}])
            # Closest clade no TR (n > 0.8, n>3.2, ...)
            attributes_divergence.append([{'Ortholog_nD': ['max',1]}])
            attributes_divergence.append([{'Ortholog_nD': ['max',0.8]}])
            # Closest clade Difference number repeat units         
            attributes_divergence.append([{'dn': ['min', 4]}])
            attributes_divergence.append([{'Ortholog_dn': ['min', 4]}])
            attributes_divergence.append([{'Ortholog_dn': ['min', 5]}])
            # Closest clade Some divergence
            attributes_divergence.append([{'dn': ['min', 4], 'd(max_n,parsimony)' : ['min', 4]}])
            # Closest clade not perfect
            attributes_divergence.append([{'dn': ['min', 1]}, {'Kendall': ['less', 1]}, {'d(max_n,nBC)': ['min', 1]}])
            for i in attributes_divergence:
                clade_distance = most_extreme_TR(paralog_Data = iParalog_Data, species_main = species_main, lConstraint = i, distant = False)
                iMerged_Data[str(i)] = clade_distance
    
    with open(result_file, 'wb') as result_handle:
        pickle.dump(all_TR, result_handle) 


def valid(constraint, data):
    
    for constraint_name, constraint_value in constraint.items():
        try:
            if constraint_name == 'd(max_n,nBC)':
                my_Value = max(data['n0'],data['n1']) - data['nBC']
            elif constraint_name == 'd(min_n,nBC)':
                my_Value = max(data['n0'],data['n1']) - data['nBC'] 
            elif constraint_name == 'd(min_n,parsimony)':
                my_Value = min(data['n0'],data['n1']) - data['min_split']
            elif constraint_name == 'd(max_n,parsimony)':
                my_Value = max(data['n0'],data['n1']) - data['min_split']
            elif constraint_name == 'Ortholog_dn':
                my_Value = abs(data['Ortholog_n'] - data['n'])
            elif constraint_name == 'Ortholog':
                return True 
            else:
                my_Value = data[constraint_name]
        except:
            #print(data)
            #print(constraint_name)
            return False
        
        if my_Value == 'NA':
            return False
        if constraint_value[0] == 'min':
            if my_Value < constraint_value[1]:
                return False
        elif constraint_value[0] == 'less':
            if my_Value <= constraint_value[1]:
                return False
        else:
            if my_Value > constraint_value[1]:
                return False
    
    return True

def most_extreme_TR(paralog_Data, species_main = 'Homo_sapiens', lConstraint = [{'nBC':0.35, 'min_Kendall':0.35}], distant = True):
        
    ''' Find most distant ortholog that confirms with lConstraint '''
    if distant:
        distance_extreme = -1
    else:
        # Warning! Works only if the maximum number of clades is 18.
        distance_extreme = 19
    
    for iO, data in paralog_Data.items():
  
        iO_clade_distance = ensembl_species_info.translate("ensembl_ID", "clade_distance",iO)
        # Check if the <clade_distance> of iP, iO is higher than the already established max distance.
        if (distant and iO_clade_distance <= distance_extreme) or (not distant and iO_clade_distance >= distance_extreme) :
            continue
        
        # Are any of the constraints in <lConstraint> met?
        if any(valid(constraint, data) for constraint in lConstraint):
            distance_extreme = iO_clade_distance
    
    return distance_extreme
    

def iterate_tandem_repeat_data(tandem_repeat_dir):
    
    lFile = [iF for iF in os.listdir(tandem_repeat_dir) if iF.endswith(".pickle")]
    lTR_ID = [iFile.split(".")[0] for iFile in lFile] 
    
    for iCount, iTR_ID in enumerate(lTR_ID):
        print("{}:{} {}".format(iCount,len(lTR_ID),iTR_ID))
        my_TR_data = gene_tree_info.load_pickles(iTR_ID, tandem_repeat_dir)
        yield iTR_ID, my_TR_data

############################## WRITE RESULT DATA ##############################

def write_tandem_repeat_data_per_repeat(data_file, result_file):
    
    with open(data_file,'rb') as dfh:
        data = pickle.load(dfh)   
    
    attributes_Human = ['lHMM', 'lMean', 'n', 'nD', 'detection_ID','nExon', 'nTR_in_Exon_Max','Ensembl_Protein_ID','divergence']
    attributes_Extreme = list(dTranslate.keys())
    
    with open(result_file, 'w') as result_handle:
        result_handle.write("\t".join(["TR_ID"]+attributes_Human+attributes_Extreme) + "\n")
        for my_TR_ID, merged_data_TR in data['merged_data'].items():
            print(my_TR_ID)
            for ensembl_ID, p_merged_data in merged_data_TR.items():
                for o_all_data in data['all_data'][my_TR_ID][ensembl_ID].values():
                    pAttributes_Human = [o_all_data[i] for i in attributes_Human]
                    break
                
                pAttributes_Extreme = [p_merged_data[dTranslate[i]] for i in attributes_Extreme]
                pAttributes =  [str(iA) for iA in pAttributes_Human + pAttributes_Extreme]
                result_handle.write("{}\t".format(str(my_TR_ID)))
                result_handle.write("\t".join(pAttributes) + "\n")


################# NOT CONSERVED TANDEM REPEATS OUTPUT #################

def non_conserved_attributes(tandem_repeat_dir, clade_distance = 0, result_file_dir_name = 'homininae_non_perfect', result_file_name = 'none_perfect_attributes.tsv'):
    
    if result_file_dir_name:
        result_dir = os.path.join(tandem_repeat_dir,result_file_dir_name)
    else:
        result_dir = tandem_repeat_dir
    
    attributes = ['lHMM','lMean', 'n', 'nD', 'detection_ID','exon', 'Ensembl_Protein_ID', 'Ensembl_Protein_ID_Ortholog', 'n1', 'dn', 'min_split', 'nC','nCB','p_nCB','dnBC_nC','Kendall', 'p_Kendall']
    short_attributes = ['lHMM','lMean', 'n', 'nD', 'detection_ID','exon', 'Ensembl_Protein_ID']
    
    dAttribute_Name_Conversion = {'exon': 'nExon\tnTR_in_Exon_Max'}
    attributes_names = ['closest_clade_name', 'closest_clade_distance'] + [dAttribute_Name_Conversion[i] if i in dAttribute_Name_Conversion else i for i in attributes]
    
    result_file = os.path.join(result_dir,result_file_name)
    
    with open(result_file, 'w') as result_handle:
        result_handle.write("TR_ID\t{}\n".format("\t".join(attributes_names)))
        for my_TR_ID, my_TR_data in iterate_tandem_repeat_data(tandem_repeat_dir):
            
            for ensembl_ID in my_TR_data['homologous_pairs'].keys():
                
                # Look only at TRs where the gene is present at <clade_distance>, but not the TR.
                conservation_data = my_TR_data['conservation']['paralog_wise'][ensembl_ID]
                if not (conservation_data['ortholog_present_oldest'][1] >= clade_distance and conservation_data['most_distant_clade']["[{'dn': 0, 'Kendall': 1, 'max_d(max_n,nBC)': 0}]"][1] < clade_distance):
                    continue
                
                print(ensembl_ID)
                # Which is the <closest_clade> to <clade_distance> in which there is a gene ortholog?
                closest_clade = (None,100)
                lOrtholog_subset = []
                for iO, iClade in zip(my_TR_data['orthology'][ensembl_ID]['lOrtholog'],my_TR_data['orthology'][ensembl_ID]['lClade']):
                    if iClade[1] < closest_clade[1] and iClade[1] >= clade_distance:
                        closest_clade = iClade
                    if iClade[1] == clade_distance:
                        lOrtholog_subset.append(iO)
                closest_clade = list(closest_clade)
                
                # In case <closest_clade> == clade_distance, look at all orthologs in that distance
                if len(lOrtholog_subset) > 0:
                    for iO in lOrtholog_subset:
                        
                        # Make files to create the tree
                        if iO in my_TR_data['homologous_pairs'][ensembl_ID]:
                            ### UNCOMMENT IF SCRIPTREE FILES NEEDED!
                            tree_io.create_scriptree_files(result_file_stump = os.path.join(result_dir, "{}_{}_{}".format(my_TR_ID,ensembl_ID,iO)), tree = my_TR_data['homologous_pairs'][ensembl_ID][iO], annotations = None)
                        
                            # Save the following features:
                            # n0, n1, nC, nBC, Kendall, gene_present_oldest,
                            iAttributes = closest_clade + get_attributes(my_TR_data, ensembl_ID, attributes, ortholog = iO)
                        else:
                            iAttributes = closest_clade + get_attributes(my_TR_data, ensembl_ID, short_attributes)
                            iAttributes = iAttributes + (len(attributes_names)-len(iAttributes) + 1)*['NA']
                        iAttributes = [str(iA) for iA in iAttributes]
                        result_handle.write("{}\t".format(str(my_TR_ID)))
                        result_handle.write("\t".join(iAttributes) + "\n")
                else:
                    iAttributes = closest_clade + get_attributes(my_TR_data, ensembl_ID, short_attributes)
                    iAttributes = iAttributes + (len(attributes_names)-len(iAttributes) + 1)*['NA']
                    iAttributes = [str(iA) for iA in iAttributes]
                    result_handle.write("{}\t".format(str(my_TR_ID)))
                    result_handle.write("\t".join(iAttributes) + "\n")
                    
                    
################# SPECIES WISE OUTPUT #################
    
def get_species_wise_info(all_TR, species = 'Pan_troglodytes', constraint = {'lD': 15, 'nD': 3.5}):
    
    ''' deprecated ?''' 
    
    species_data = []
    for my_TR_ID, my_TR_data in all_TR.items():
        for my_paralog_data in my_TR_data.values():
            for iO, my_ortholog_data in my_paralog_data.items():
                if species != ensembl_species_info.translate("ensembl_ID", "official",iO):
                    continue
                if any(float(my_ortholog_data[i]) < j for i,j in constraint.items()):
                    break
                my_ortholog_data["TR_ID"] = my_TR_ID
                my_ortholog_data['lD'] = int(my_ortholog_data['lD'])
                species_data.append(my_ortholog_data)
    return species_data           

    
################################## STATISTICS ##################################   


def get_IDs(data_file):
    ''' Get all used Ensembl IDs ordered by species '''
    
    with open(data_file,'rb') as dfh:
        all_TR = pickle.load(dfh)
    
    # Save lists of all Ensembl_IDs that comply with the criterium
    result = defaultdict(list)
    for TR_ID, TR_data in all_TR['merged_data'].items():
        for ensembl_ID, p_data in TR_data.items():
            result[ensembl_species_info.translate('ensembl_ID','english', ensembl_ID)].append(ensembl_ID)
            for o_ensembl_ID in all_TR['all_data'][TR_ID][ensembl_ID].keys():
                result[ensembl_species_info.translate('ensembl_ID','english', o_ensembl_ID)].append(o_ensembl_ID)
    return result
    
def get_subset(data_file, result_dir, criterium = ['furthest_TR_p', 8, 'min']):
    
    ''' Produce lists of EnsemblIDs and PFAM_IDs for TRs with certain characteristics, 
        such as
        - Conserved in Eutheria
        - Diverged in Eutheria
        - Save these lists, ordered by PFAM_ID
        
        criterium in ['furthest_TR_p','furthest_TR_OG_d1_K1', 'closest_TR_completely_diverged', 'closest_TR_strongly_diverged_1', 'closest_min_dn_4']
        
        - Calculate Which PFAM_IDs are particularly prominent
        - Save second list with PFAM_IDS ordered by pValue (hypergeometrical distribution)
        '''
    with open(data_file,'rb') as dfh:
        all_TR = pickle.load(dfh)
        
    with open(os.path.join(DATAROOT2, "PDB_Ensembl_Protein_ID.tsv"), "r") as fh:
        tmp = fh.readlines()
        tmp = [i.rstrip().split("\t") for i in tmp]
        dEnsembl_Protein_ID_PDB = {i[1]:i[0] for i in tmp if i[0] != ""}
    
    clade_name = ensembl_species_info.translate('clade_distance','clade_name', criterium[1])
    result_stump = os.path.join(result_dir,"{}_{}_".format(clade_name, str(criterium[0])))
    background_stump = os.path.join(result_dir,"{}_Background_".format(str(criterium[1])))
    
    # Save lists of all Ensembl_IDs that comply with the criterium
    lEnsembl_ID = []
    lEnsembl_ID_All = []
    lPFAM = []
    lPFAM_All = []
    for TR_ID, TR_data in all_TR['merged_data'].items():
        for ensembl_ID, p_data in TR_data.items():
            for o_data in all_TR['all_data'][TR_ID][ensembl_ID].values():
                pfam = o_data['detection_ID']
                lHMM = o_data['lHMM']
                nD = o_data['nD']
                break
            if nD < 3.5 or lHMM < 15:
                continue
            lPFAM_All.append(pfam)
            lEnsembl_ID_All.append(ensembl_ID)
            if (criterium[2] == 'min' and p_data[dTranslate[criterium[0]]] >= criterium[1]) or (criterium[2] == 'max' and p_data[dTranslate[criterium[0]]] <= criterium[1]):
                lEnsembl_ID.append(ensembl_ID)
                lPFAM.append(pfam)
    
    with open(result_stump + "Ensembl_PFAM_PDB.txt",'w') as rfh:
        for ensembl_ID,pfam in zip(lEnsembl_ID, lPFAM):
            rfh.write("{}\t{}\t{}\n".format(ensembl_ID, pfam, dEnsembl_Protein_ID_PDB[ensembl_ID] if ensembl_ID in dEnsembl_Protein_ID_PDB else ""))
    with open(background_stump + "Ensembl_PFAM_PDB.txt",'w') as rfh:
        for ensembl_ID,pfam in zip(lEnsembl_ID_All, lPFAM_All):
            rfh.write("{}\t{}\t{}\n".format(ensembl_ID, pfam, dEnsembl_Protein_ID_PDB[ensembl_ID] if ensembl_ID in dEnsembl_Protein_ID_PDB else ""))
    with open(result_stump + "Ensembl.txt",'w') as rfh:
        rfh.write("\n".join(lEnsembl_ID))
    with open(result_stump + "INVERSE_Ensembl.txt",'w') as rfh:
        rfh.write("\n".join([i for i in lEnsembl_ID_All if i not in lEnsembl_ID]))
    with open(background_stump + "Ensembl.txt",'w') as rfh:
        rfh.write("\n".join(lEnsembl_ID_All))
    
    # Do enrichment analysis of PFAM terms (on lPFAM, lPFAM_All)
    dPFAM = defaultdict(int)
    for i in lPFAM:
        dPFAM[i] += 1
    
    dPFAM_All = defaultdict(int)
    for i in lPFAM_All:
        dPFAM_All[i] += 1
    
    dSignificance = defaultdict(int)
    for i in dPFAM.keys():
        dSignificance[i] = get_significance(dPFAM_All, dPFAM, i)
    
    # Ordered by p-Value
    lData = []
    with open(result_stump + "enriched.txt",'w') as rfh:
        for iPFAM, iSig in sorted(dSignificance.items(), key=lambda x: x[1]):
            rfh.write("{}\t{}\t{}\t{}\n".format(iPFAM, iSig, dPFAM[iPFAM], dPFAM_All[iPFAM]))
            pfam_name = dPFAM_Name[iPFAM] if iPFAM in dPFAM_Name else '???'
            lData.append([pfam_name, iPFAM, "{}/{}".format(dPFAM[iPFAM],dPFAM_All[iPFAM]), "${0:.1e}}}$".format(iSig).replace("e-0"," \cdot 10^{-").replace("e-"," \cdot 10^{-")])
    # Produce Table 3 in Publication.
    # Description PFAM count p-Value       
    lTable_Data = [" & ".join(["Description", "PFAM ID", "Conserved", "pValue"]) + " \\\\"]
    lTable_Data += [" & ".join(iData) + "\\\\" for iData in lData]
    
    # Ordered by value
    lData = []
    dAbsolut = {iPFAM: (dPFAM[iPFAM]/dPFAM_All[iPFAM],dPFAM[iPFAM]) for iPFAM in dPFAM.keys()}
    with open(result_stump + "enriched_absolute.txt",'w') as rfh:
        for iPFAM, iSig in sorted(dAbsolut.items(), key=lambda x: x[1], reverse = True):
            rfh.write("{}\t{}\t{}\t{}\n".format(iPFAM, dSignificance[iPFAM], dPFAM[iPFAM], dPFAM_All[iPFAM]))
            pfam_name = dPFAM_Name[iPFAM] if iPFAM in dPFAM_Name else '???'
            lData.append([pfam_name, iPFAM, "{}/{}".format(dPFAM[iPFAM],dPFAM_All[iPFAM]), "${0:.1e}}}$".format(dSignificance[iPFAM]).replace("e-0"," \cdot 10^{-").replace("e-"," \cdot 10^{-")])
    # Produce Table 3 in Publication.
    # Description PFAM count p-Value       
    lTable_Data_abs = [" & ".join(["Description", "PFAM ID", "Conserved", "pValue"]) + " \\\\"]
    lTable_Data_abs += [" & ".join(iData) + "\\\\" for iData in lData]    
    
    return lTable_Data,lTable_Data_abs

def get_TR_evolution_summary_table_publication(data_file, result_dir):
    
    ''' Produce Table 2 in Publication.
    E.g. extract general, species wise and clade wise data from your TR summary file <data_file>
        '''
    
    lSpecies = ['Chimpanzee','Mouse', 'Rat', 'Xenopus', 'Fugu', 'Zebrafish', 'Fruitfly', 'Roundworm', "Yeast"]    
    lClade = [0,8,14]
    lOrder = [(True, 0),(False, 0),(True, 1),(False, 1),(False, 2),(True, 2),(False, 3),(False, 4),(False, 5),(False, 6),(False, 7),(False, 8)]
    #lSubset = ['all','PF00096', 'PF00560', 'PF00400', 'PF00023', 'PF00028', 'PF07679', 'PF07645', 'PF01344', 'de novo']
    lSubset = ['all','PF00096', 'PF00560', 'PF00400', 'PF00023', 'PF00028', 'PF07679', 'de novo']
    dConstraint = {'perfectly_conserved': {'dn': ['max', 0], 'Kendall': ['min', 1], 'd(max_n,nBC)': ['max', 0]}, 'completely_diverged': {'n1': ['min', 4], 'min_split': ['max', 1]}}
    lConstraint = ['perfectly_conserved','completely_diverged']
     
    with open(data_file,'rb') as dfh:
        all_TR = pickle.load(dfh)
    
    # Get general information
    lData = ['all','mean(l)','mean(n)','std(n)','mean(div)']
    dData = {'mean(l)':('lMean',np.mean), 'mean(n)':('nD',np.mean),'std(n)':('nD',np.std), 'mean(div)':('divergence',np.mean)}
    
    data_results = {iData: {iSubset: [] for iSubset in lSubset} for iData in lData}
    data_results['all'] = {iSubset: [] for iSubset in lSubset}
    lDivergence = {iSubset: [] for iSubset in lSubset}
    for TR_ID, TR_data in all_TR['all_data'].items():
        for ensembl_ID, p_data in TR_data.items():
            for o_data in p_data.values():
                pfam = o_data['detection_ID']
                lHMM = o_data['lHMM']
                nD = o_data['nD']
                break
            if nD < 3.5 or lHMM < 15:
                continue
            data_results['all']['all'].append(1)
            if pfam in lSubset:
                data_results['all'][pfam].append(1)
            elif pfam in ['hhrepid', 'xstream', 't-reks']:
                data_results['all']['de novo'].append(1)
            for iData,iData_info in dData.items():
                data_results[iData]['all'].append(o_data[iData_info[0]])
                if pfam in lSubset:
                    data_results[iData][pfam].append(o_data[iData_info[0]])
                elif pfam in ['hhrepid', 'xstream', 't-reks']:
                    data_results[iData]['de novo'].append(o_data[iData_info[0]])
    
    lTable_Data = []
    for iData in lData:
        tmp = [iData]
        for iSubset in lSubset:
            if iData == 'all':
                tmp2 = "{}".format(np.sum(data_results[iData][iSubset]))
            elif iData == 'mean(div)':
                tmp2 = "{:.2f}".format(float(dData[iData][1](data_results[iData][iSubset])))
            else:
                tmp2 = "{:.1f}".format(float(dData[iData][1](data_results[iData][iSubset])))
            if iSubset == 'de novo':
                tmp.append("\multicolumn{3}{|c}{" + tmp2 + "}")
            elif iSubset == 'all':
                tmp.append("\multicolumn{3}{c|}{" + tmp2 + "}")
            else:
                tmp.append("\multicolumn{3}{c}{" + tmp2 + "}")
        lTable_Data.append(" & ".join(tmp) + "\\\\")
    
    # Get clade-wise information.
    clade_results = {iClade: {iSubset:{'all':0, 'perfectly_conserved':{'count': 0, 'pValue':None}, 'completely_diverged':{'count': 0, 'pValue':None}} for iSubset in lSubset} for iClade in lClade}
    for TR_ID, TR_data in all_TR['merged_data'].items():
        for ensembl_ID, p_data in TR_data.items():
            for o_data in all_TR['all_data'][TR_ID][ensembl_ID].values():
                pfam = o_data['detection_ID']
                lHMM = o_data['lHMM']
                nD = o_data['nD']
                break
            if nD < 3.5 or lHMM < 15:
                continue
            for iClade in lClade:
                if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_n_4']] >= iClade:
                    clade_results[iClade]['all']['all'] += 1
                if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_p']] >= iClade:
                    clade_results[iClade]['all']['perfectly_conserved']['count'] += 1
                if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['closest_TR_diverged']] <= iClade:
                    clade_results[iClade]['all']['completely_diverged']['count'] += 1
                if pfam in lSubset:   
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_n_4']] >= iClade:
                        clade_results[iClade][pfam]['all'] += 1
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_p']] >= iClade:
                        clade_results[iClade][pfam]['perfectly_conserved']['count'] += 1
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['closest_TR_diverged']] <= iClade:
                        clade_results[iClade][pfam]['completely_diverged']['count'] += 1
                elif pfam in ['hhrepid', 'xstream', 't-reks']:  
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_n_4']] >= iClade:
                        clade_results[iClade]['de novo']['all'] += 1
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['furthest_TR_p']] >= iClade:
                        clade_results[iClade]['de novo']['perfectly_conserved']['count'] += 1
                    if all_TR['merged_data'][TR_ID][ensembl_ID][dTranslate['closest_TR_diverged']] <= iClade:
                        clade_results[iClade]['de novo']['completely_diverged']['count'] += 1                    
    
    # Create Table rows for clade data.
    lTable_Clade = []
    for iClade in lClade:
        tmp = [ensembl_species_info.translate("clade_distance", "clade_name",iClade)]
        for iSubset in lSubset:
            lConstrained_Count = []
            for iConstraint in lConstraint:
                lConstrained_Count.append(str(clade_results[iClade][iSubset][iConstraint]['count']))
            count = str(clade_results[iClade][iSubset]['all'])
            tmp.append("{} &({},&{})".format(count, *lConstrained_Count)) 
        lTable_Clade.append("\\hline \\rowcolor{lightgray} " + " & ".join(tmp) + "\\\\")
        
    # Get species-wise information.
    species_results = {iSpecies: {iSubset:{'all':0, 'perfectly_conserved':{'count': 0, 'pValue':None}, 'completely_diverged':{'count': 0, 'pValue':None}} for iSubset in lSubset} for iSpecies in lSpecies}
    
    dN = {iSpecies: {iSubset: defaultdict(int) for iSubset in lSubset} for iSpecies in lSpecies}
    for TR_ID, TR_data in all_TR['all_data'].items():
        for ensembl_ID, p_data in TR_data.items():
            # Find all orthologs to <ensembl_ID> in <lSpecies> & See how many are conserved
            for o_ensembl_ID, o_data in p_data.items():
                lHMM = o_data['lHMM']
                nD = o_data['nD']
                if nD < 3.5 or lHMM < 15:
                    break
                ortholog_nD = o_data['Ortholog_nD']   
                iSpecies = ensembl_species_info.translate('ensembl_ID','english', o_ensembl_ID)
                if ortholog_nD < 3.2 or iSpecies not in lSpecies:
                    continue
                pfam = o_data['detection_ID']
                tN = (min(o_data['n0'],o_data['n1']), max(o_data['n0'],o_data['n1']))
                if tN[0] == 'NA' or tN[1] == 'NA':
                    print("NA: {} {} {}".format(TR_ID, ensembl_ID, o_ensembl_ID))
                    continue
                dN[iSpecies]['all'][tN] += 1
                species_results[iSpecies]['all']['all'] += 1
                if pfam in lSubset:               
                    species_results[iSpecies][pfam]['all'] += 1
                    dN[iSpecies][pfam][tN] += 1
                elif pfam in ['hhrepid', 'xstream', 't-reks']:
                    species_results[iSpecies]['de novo']['all'] += 1
                    dN[iSpecies]['de novo'][tN] += 1
                for constraint_name, constraint in dConstraint.items():
                    if valid(constraint, o_data):
                        species_results[iSpecies]['all'][constraint_name]['count'] += 1
                        if pfam in lSubset:
                            species_results[iSpecies][pfam][constraint_name]['count'] += 1
                        if pfam in ['hhrepid', 'xstream', 't-reks']:
                            species_results[iSpecies]['de novo'][constraint_name]['count'] += 1
    
    # Calculate pValue of the species_wise results
    for iSpecies, species_data in species_results.items():
        for iSubset, subset_data in species_data.items():
            p_combinations = dN[iSpecies][iSubset]
            # Conserved
            nObserved = subset_data['perfectly_conserved']['count']
            try:
                nExpected, pValue = get_expected_conserved(p_combinations, nObserved)
            except:
                print("{} {}".format(iSpecies, iSubset))
                print(nObserved)
                print(p_combinations)
            subset_data['perfectly_conserved']['pValue'] = pValue
            # Diverged
            nObserved = subset_data['completely_diverged']['count']
            try:
                nExpected, pValue = get_expected_completely_diverged(p_combinations, nObserved)
            except:
                print("{} {}".format(iSpecies, iSubset))
                print(nObserved)
                print(p_combinations)
            subset_data['completely_diverged']['pValue'] = pValue
    
    # Create Table rows for species data.
    lTable_Species = []
    for iSpecies in lSpecies:
        tmp = [iSpecies]
        for iSubset in lSubset:
            lConstrained_Count = []
            for iConstraint in lConstraint:
                count = str(species_results[iSpecies][iSubset][iConstraint]['count'])
                if species_results[iSpecies][iSubset][iConstraint]['pValue'] <= 0.01:
                    count = "\textcolor{Clover}{" + count + "}"
                elif species_results[iSpecies][iSubset][iConstraint]['pValue'] >= 0.99:
                    count = "\textcolor{Maraschino}{" + count + "}"
                lConstrained_Count.append(count)
            count = str(species_results[iSpecies][iSubset]['all'])
            tmp.append("{} &({},&{})".format(count, *lConstrained_Count))
        lTable_Species.append(" & ".join(tmp) + "\\\\")
    
    allTable = lTable_Data + [lTable_Clade[iOrder[1]] if iOrder[0] else lTable_Species[iOrder[1]] for iOrder in lOrder]
    return(allTable)   
        

def get_significance(all, partition, key):
    
    nAll = sum(all.values())
    nPartition = sum(partition.values())
    nPFAM_ALL = all[key]
    nPFAM_Partition = partition[key]
    
    if nPFAM_Partition == 0:
        return 1
    
    return hypergeom.sf(nPFAM_Partition-1,nAll,nPFAM_ALL,nPartition)    
    
    
def get_expected_conserved(p_combinations, nObserved):
    
    ''' Calculate how many perfectly conserved TRs you expect within <p_combinations> and  <nObserved>.
        
        p_combinations = {(nA,nB): count}
        e.g.
        p_combinations = {(2,2): 4, (2,3):2, ...}
        
        <nObserved> is the number of actually observed conserved TRs.
    '''
    
    p_possibly_perfect = {i[0]:j for i,j in p_combinations.items() if i[0]==i[1]}
    
    if len(p_possibly_perfect) == 0:
        return 0, 1
    
    # Get probability for each value of n.
    p_distributions = [[binom.pmf(iCount,count,repeat_tree_info.probability_perfect_conservation(n)) for iCount in range(count+1)] if n < 40 else [1]+[0]*(count) for n,count in p_possibly_perfect.items()]
    #return p_distributions
    
    # Calculate the convolved distribution.
    p_convolved = p_distributions[0]
    for i in p_distributions[1:]:
        p_convolved = sc.convolve(p_convolved, i)
    
    ##return p_convolved, p_distributions, p_possibly_perfect
    # Calculate the number of expected perfect cases.
    nExpected = sum(i*j for i,j in enumerate(p_convolved))
    
    # What's the pValue?
    pValue = sum(p_convolved[nObserved:])
    
    return nExpected, pValue
    
def get_expected_completely_diverged(p_combinations, nObserved):
       
    ''' Calculate how many perfectly diverged TRs you expect within <p_combinations> and  <nObserved>.'''   
    
    # Get probability for each value of n.
    p_distributions = [[binom.pmf(iCount,count,repeat_tree_info.probability_completely_diverged(n[0],n[1])) for iCount in range(count+1)] if repeat_tree_info.probability_completely_diverged(n[0],n[1]) > 0 else [1]+[0]*(count) for n,count in p_combinations.items()]
        
    # Calculate the convolved distribution.
    p_convolved = p_distributions[0]
    for i in p_distributions[1:]:
        p_convolved = sc.convolve(p_convolved, i)
    
    # Calculate the number of expected perfect cases.
    nExpected = sum(i*j for i,j in enumerate(p_convolved))
    
    # What's the pValue?
    pValue = sum(p_convolved[nObserved:])
    
    return nExpected, pValue
   
################# CREATE MULTITREE #################
   
def create_multitree(tandem_repeat_dir, gene_tree_dir, tandem_repeat_ID, ensembl_ID, result_dir, lSpecies = ['Mus_musculus']):
   
    my_TR_data, my_gene_tree = gene_tree_info.load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)
    lID = [ensembl_ID]
    lID += [iO for iO in my_TR_data['orthology'][ensembl_ID]['lOrtholog'] if ensembl_species_info.translate("ensembl_ID", "official",iO) in lSpecies]
    aligned_msas, newick_tree = create_tree(tandem_repeat_ID, my_TR_data, my_gene_tree, lID, result_dir, key = None, create_draw_tree_files = True)
    
    
################# EXON OUTPUT #################

def write_exon_data(tandem_repeat_dir, gene_tree_dir, result_file_name = 'exon.ex'):
    
    ''' Visualise the exon data for all TR pickles in <tandem_repeat_dir>'''
    
    result_file = os.path.join(tandem_repeat_dir,result_file_name)
    result_file2 = os.path.join(tandem_repeat_dir,result_file_name+"2")
    lFile = [iF for iF in os.listdir(tandem_repeat_dir) if iF.endswith(".pickle")]
    lTR_ID = [iFile.split(".")[0] for iFile in lFile]      
    
    
    with open(result_file2, 'w') as result_handle2:
        with open(result_file, 'w') as result_handle:
            for iTR_ID in lTR_ID:
                print(iTR_ID)
                my_TR_data, my_gene_tree = gene_tree_info.load_pickles(iTR_ID, tandem_repeat_dir, gene_tree_dir)
                for ensembl_ID, data in my_TR_data['homologous_pairs'].items():
                    exon_structure = my_gene_tree['leafs'][ensembl_ID]['exon']
                    my_TR = my_TR_data['homologs'][ensembl_ID]['TR']['HMM']
                    if type(my_TR) == dict:
                        continue
                    exon = repeat_exon.print_exon_structure(my_TR, exon_structure)            
                    result_handle.write(">{} {} \n{}\n".format(iTR_ID, ensembl_ID, exon))
                    result_handle2.write(">{} {} \n{}\n".format(iTR_ID, ensembl_ID, exon.replace("_","")))

if __name__ == "__main__":       
    
    print('Done')
    
    

        
