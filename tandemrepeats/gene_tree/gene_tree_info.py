# Copyright (C) 2012-2013 by Elke Schaper (elke@inf.ethz.ch) 

"""
    Compile pickled/XML Gene Tree and Tandem Repeat Data from 
    - Ensembl Gene Trees
    - Ensembl Sequence Information
    - Tandem Repeat Detection
    
    Conduct Homolog Tandem Repeat Detection on Gene Tree data
    
""" 

import itertools,logging,os,pickle,re
from collections import defaultdict
from Bio import SeqIO
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.propagate = False
#logging.basicConfig(level=logging.INFO)

from . import align, ensembl_IO, tree_io
from repeat.hmm import hmm, viterbi
from repeat.paths import *
from ..repeat import repeat_info, repeat_realign, repeat_io

    
def assemble_gene_tree_data(gene_tree, gene_tree_ID, result_file):

    ''' Add (mostly sequence) data to <gene_tree> dict.'''

    ## UNCOMMENT! (takes loads of time)
    record_dict_compara_aa = SeqIO.index(os.path.join(DATAROOT2,'Compara.69.protein.aa.fasta'), "fasta")
    #record_dict_compara_cds = SeqIO.index(os.path.join(DATAROOT2,'Compara.69.protein.cds.fasta'), "fasta")

    logging.info("Gene_tree: {0}".format(gene_tree_ID))
    # gene_tree is a dictionary with a dictionary of all leaf sequences <leafs> and the tree in nhx format: <tree_nhx>
    # Save <gene_tree_ID>, <ancestral_clade> (unknown so far) to gene_tree:
    gene_tree['gene_tree_ID'] = gene_tree_ID
    gene_tree['ancestral_clade'] = None
    
    iLeaf = 0
    for ensembl_protein_ID,leaf in gene_tree['leafs'].items():
        iLeaf += 1
        # leaf contains the values for ensembl_gene_ID, species
        logging.debug("ensembl_protein_ID in gene_tree: {0}".format(ensembl_protein_ID))
        logging.debug("Leaf in gene_tree: {0}".format(str(leaf)))
        species = leaf['species']
        logging.info("Species in gene_tree: {0}".format(species))
        ensembl_gene_ID = leaf['gene_ID']
        logging.debug("ensembl_gene_ID in gene_tree: {0}".format(ensembl_gene_ID))
        
        leaf['sequences'] = {}
        
        ## Save gene tree alignment protein and coding sequence UNCOMMENT
        leaf['sequences']['aa_gene_tree_alignment'] = record_dict_compara_aa[ensembl_protein_ID].seq
        #leaf['sequences']['cds_gene_tree_alignment'] = record_dict_compara_cds[ensembl_protein_ID].seq  
        
        # Save Ensembl Transcript ID, GO Term, Transcript count, Status(Gene), Status(transcript), GC content from '<species>_Xrefs.txt'
        ## GO: There may be several GO Terms per entry. At current, you are only reading the first line = the first entry
        logging.debug("xrefs File: {0}".format(os.path.join(DATAROOT2,'_'.join([species, 'xrefs.txt']))))
        with open(os.path.join(DATAROOT2,'_'.join([species, 'xrefs.txt'])), 'r') as xref_filehandle:
            xrefs = ensembl_IO.get_data(xref_filehandle, ("Ensembl Protein ID",ensembl_protein_ID))
        # Later: Transfer further Xrefs to the dictionary leaf.
        logging.debug("Xrefs: {0}".format(str(xrefs)))
        leaf['ensembl_transcript_ID'] = xrefs['Ensembl Transcript ID']
        leaf['gene_start_bp'] = xrefs['Gene Start (bp)'] 
        leaf['gene_end_bp'] = xrefs['Gene End (bp)'] 
        ensembl_transcript_ID = leaf['ensembl_transcript_ID']
        logging.debug("ensembl_transcript_ID: {0}".format(ensembl_transcript_ID))
        
        # Save coding sequence from '<species>_cds.faa' via <ensembl_transcript_ID>
        logging.debug("CDS File: {0}".format(os.path.join(DATAROOT2,'_'.join([species, 'cds.fasta']))))
        record_dict_cds = SeqIO.index(os.path.join(DATAROOT2,'_'.join([species, 'cds.fasta'])), "fasta")
        try:
            leaf['sequences']['cds'] = record_dict_cds[ensembl_protein_ID]
        except:
            logging.error("Could not detect coding sequence (cds) for {0}".format(ensembl_protein_ID)) 
        

        # Save protein sequence from '<species>_aa.faa' via <ensembl_transcript_ID>
        logging.debug("aa File: {0}".format(os.path.join(DATAROOT2,'_'.join([species, 'aa.fasta']))))
        record_dict_aa = SeqIO.index(os.path.join(DATAROOT2,'_'.join([species, 'aa.fasta'])), "fasta")
        try:
            leaf['sequences']['aa'] = record_dict_aa[ensembl_protein_ID]
        except:
            logging.error("Could not detect protein sequence (aa) for {0}".format(ensembl_protein_ID)) 

        # Save unspliced gene sequence from '<species>_unspliced_gene.faa' via <ensembl_transcript_ID>
        logging.debug("Unspliced Gene File: {0}".format(os.path.join(DATAROOT2,'_'.join([species, 'unspliced_gene.fasta']))))
        record_dict_unspliced_gene = SeqIO.index(os.path.join(DATAROOT2,'_'.join([species, 'unspliced_gene.fasta'])), "fasta")
        try:
            leaf['sequences']['unspliced_gene'] = record_dict_unspliced_gene[ensembl_gene_ID]
        except:
            logging.error("Could not detect unspliced_gene sequence for {0}".format(ensembl_gene_ID)) 
        
        
        # FOR EACH Exon:
        ## Save Exon Chr Start (bp), Exon Chr End (bp), CDS Start, CDS End from '<species>_exon_structure.txt'
        ## If you want to rename the exon dict keys, you'll have to change leaf['exons']
        logging.debug("Exon File: {0}".format(os.path.join(DATAROOT2,'_'.join([species, 'exon.txt']))))
        with open(os.path.join(DATAROOT2,'_'.join([species, 'exon.txt'])), 'r') as exon_filehandle:
            leaf['exon'] = ensembl_IO.get_data(exon_filehandle, ("Ensembl Protein ID",ensembl_protein_ID), all = True)
         
        logging.debug("exons: {0}".format(leaf['exon'])) 
        leaf['tandem repeats'] = defaultdict(list)
 
                
## Once you decide to look at nucleic TRs also, decide at what stage to translate cds into aa coordinates, or vice versa.
## (Maybe you should define your coordinate system in the beginning, and have translating functions.) 
#         # Run Nucleic TRDs on Unsliced Gene Sequence
#         # Define unsliced_gene_sequence_records as a list of Bio.Seq objects
#         unsliced_gene_sequence_records = [leaf['sequences']['unspliced_gene']]
#         nucleic_tandem_repeats = brutus_worker.runTRD(unsliced_gene_sequence_records, sequence_type = 'AA', default = True, classifiers = ['phylo', 'phylo_gap01', 'phylo_gap001'], l = False, n = False, calculate_pValues = True, calculate_scores = True, calc_divergence = True)       
#         for iTRD,sequence in nucleic_tandem_repeats.items():
#             leaf['tandem repeats']['cds'].append(sequence[0]) # As there is only one sequence, the first member in the sequences list contains all interesting information
#             for iR in leaf['tandem repeats']['cds'][-1]:
#                 # Save l, n, MSA, TRD, scores, sequence_type, position in sequence of given type
#                 iR.TRD = iTRD  
#                 # Save ID
#                 iR.ID = count_tandem_repeat
#                 count_tandem_repeat += 1
#                 # Save start and end position in alignment     
#                 iR.position = repeat_info.calculate_position_in_alignment(iR.begin, iR.sequence_length,leaf['sequences']['cds_gene_tree_alignment']) 
                    
    logging.debug("Result file: {0}".format(result_file))
    with open(result_file, 'wb') as result_filehandle:
        pickle.dump(gene_tree, result_filehandle)
    
    return gene_tree
    
    
def assemble_exon_data(gene_tree_file):
    
    ''' add exon data to gene tree pickle '''
    
    with open(gene_tree_file, 'rb') as gene_tree_filehandle:
        gene_tree = pickle.load(gene_tree_filehandle)
    
    for ensembl_protein_ID,leaf in gene_tree['leafs'].items():
        # FOR EACH Exon:
        ## Save Exon Chr Start (bp), Exon Chr End (bp), CDS Start, CDS End from '<species>_exon_structure.txt'
        ## If you want to rename the exon dict keys, you'll have to change leaf['exons']
        logging.debug("Exon File: {0}".format(os.path.join(DATAROOT2,'_'.join([leaf['species'], 'exon.txt']))))
        with open(os.path.join(DATAROOT2,'_'.join([leaf['species'], 'exon.txt'])), 'r') as exon_filehandle:
            exon_data = ensembl_IO.get_data(exon_filehandle, ("Ensembl Protein ID",ensembl_protein_ID), all = True) 
            for iExon in exon_data:
                del iExon['Ensembl Protein ID']
            leaf['exon'] = exon_data
    
    with open(gene_tree_file, 'wb') as gene_tree_filehandle:
        pickle.dump(gene_tree, gene_tree_filehandle)



def predict_tandem_repeats_on_gene_tree_data(gene_tree_ID,count_tandem_repeat):
    
    ''' Predict Tandem repeats on all sequences in a <gene tree> dictionary. '''
    
    result_file = os.path.join(RESULTROOT,"gene_tree",str(gene_tree_ID) + ".pickle")
    logging.debug("Result file: {0}".format(result_file))
    with open(result_file, 'rb') as result_filehandle:
        my_gene_tree = pickle.load(result_filehandle)
        
    for ensembl_protein_ID,leaf in my_gene_tree['leafs'].items():    

        # Define protein_sequence_records as a list of Bio.Seq objects
        protein_sequence_records = [leaf['sequences']['aa']]
        # Run Protein TRDs on Protein Sequence
        classifiers = ['phylo', 'phylo_gap01', 'phylo_gap001']
        protein_tandem_repeats = predict_tandem_repeats(protein_sequence_records, sequence_type = 'AA', default = True, classifiers = classifiers, l = False, n = False, calculate_pValues = True, calculate_scores = True)
        for iTRD,sequence in protein_tandem_repeats.items():
            # As there is only one sequence, the first member in the sequences list contains all interesting information
            for iR in sequence[0]:
                # Consider only tandem repeats that have a repeat unit predicted to be at least one character long.
                if IR.lD > 0:
                    # Save l, n, MSA, TRD, scores, sequence_type, position in sequence of given type
                    iR.TRD = iTRD  
                    # Save ID
                    iR.ID = count_tandem_repeat
                    count_tandem_repeat += 1
                    # Save start and end position in alignment
                    iR.position = {}
                    iR.position['aa_alignment'] = repeat_info.calculate_position_in_alignment(iR.begin, iR.sequence_length,leaf['sequences']['aa_gene_tree_alignment'])  
                    iR.position['aa_protein_sequence'] = {'begin': iR.begin-1, 'end': iR.begin + iR.sequence_length - 2}
                    leaf['tandem repeats']['aa'].append(iR)        

    result_file = os.path.join(RESULTROOT,"gene_tree",my_gene_tree['gene_tree_ID'] + ".pickle")
    with open(result_file, 'wb') as result_filehandle:
        pickle.dump(my_gene_tree, result_filehandle)
    
    return my_gene_tree, count_tandem_repeat
    
def predict_pfamA_tandem_repeats_on_gene_tree_data(gene_tree_ID, lPfam_ID, data_file, result_file):


    ''' 
    1. Load the pfam HMM
    2. Create internal HMM from pfam HMM (keep emission probabilities from file)
    3. Run the HMM on all human sequence in gene_tree_ID
    4. Save results.  
    '''
    
    ### Define hmmer file
    hmm_file = os.path.join(DATAROOT2, "Pfam-A.hmm")

    ### Load HMM emission probabilities from file, adapt to local HMM model.
    
    lHMM = [hmm.HMM(hmm_file=hmm_file, accession=pfam_ID) for pfam_ID in lPfam_ID]
    
    with open(data_file, 'rb') as dfh:
        my_gene_tree = pickle.load(dfh)
        
    for ensembl_protein_ID,leaf in my_gene_tree['leafs'].items():
        if not ensembl_protein_ID.startswith("ENSP0"):
            continue
                    
        print(ensembl_protein_ID)

        leaf['pfam repeats'] = []
        
        for pfam_ID,my_hmm in zip(lPfam_ID,lHMM):
            try:
                lD = my_hmm.lD
            except:
                print("HMM {} was not build.".format(pfam_ID))
                continue
                
            #if lD > 100:
            #    continue
            print("****")
            print(pfam_ID)
            
            # Define protein_sequence_records as a list of Bio.Seq objects
            protein_sequence = str(leaf['sequences']['aa'].seq).replace("*","").replace("X","")
            # Run Protein TRDs on Protein Sequence
            my_Viterbi = viterbi.Viterbi(my_hmm, protein_sequence)        
            my_Most_likely_path = my_Viterbi.viterbi()
            unaligned_msa = viterbi.hmm_path_to_non_aligned_tandem_repeat_units(protein_sequence, my_Most_likely_path, lD)
            if unaligned_msa == None or len(unaligned_msa) == 1:
                continue 
            # Realign the repeat units
            aligned_msa = repeat_realign.realign_repeat(unaligned_msa)
            my_TR = repeat_info.Repeat(aligned_msa, calc_score = True, calc_pValue = True)
            my_TR.begin = my_Most_likely_path.count('N') + 1
            my_TR.viterbi_path = my_Most_likely_path 
        
            # Save start and end position in alignment.
            my_TR.position = {}
            my_TR.position['alignment'] = repeat_info.calculate_position_in_alignment(my_TR.begin, my_TR.sequence_length,alignment = leaf['sequences']['aa_gene_tree_alignment'])  
            my_TR.position['sequence'] = {'begin': my_TR.begin-1, 'end': my_TR.begin + my_TR.sequence_length - 1}
            my_TR.Pfam = pfam_ID
            my_TR.Pfam_l = lD
            
            if 'pfam repeats' in leaf:
                leaf['pfam repeats'].append(my_TR)
            else:
                leaf['pfam repeats'] = [my_TR]    
    
    with open(result_file, 'wb') as rfh:
        pickle.dump(my_gene_tree, rfh)



def extract_tandem_repeats_per_seq_with_regard_to_pfam_HMM(gene_tree_file, lPfam_ID, result_dir, dExisting = False):
    
    '''
    1. For each pfam HMM, save tr file
    2. Are there any more (non-overlapping) TRs? IF yes, save tr file
    
    '''
   
    coordinates_blocked=[]
    # Set TR count to zero
    count_TR = 0  
    if dExisting:
        #find coordinates_blocked
        for iCount, pickle_file in dExisting.items():
            try:
                with open(pickle_file,'rb') as fh:
                    iTR = pickle.load(fh)
                # If Pfam:
                if 'PFAM_ID' in iTR['HMM']:
                    begin_alignment = min(my_TR.position['alignment']['begin'] for my_TR in iTR['homologs'].values())
                    end_alignment = max(my_TR.position['alignment']['end'] for my_TR in iTR['homologs'].values())
                    coordinates_blocked.append((begin_alignment,end_alignment))
                    count_TR += 1
                else:
                    begin_alignment = iTR['HMM']['tandem_repeat_seed'].position['aa_alignment']['begin']
                    end_alignment = iTR['HMM']['tandem_repeat_seed'].position['aa_alignment']['end']
                    #coordinates_blocked.append((begin_alignment,end_alignment))
            except:
                print("Problem - probably file could not be openend")
                count_TR += 1
    
    # Load Gene tree.
    with open(gene_tree_file, 'rb') as gh:
        my_gene_tree = pickle.load(gh)
    
    # Set parameters
    c = {'alpha' : 0.01, 'classifier' : 'phylo_gap001_ignore_trailing_gaps_and_coherent_deletions', 'lD' : 10, 'd': 0.8, 'nD':2.5} 

    # Find homo sapiens within species <lParalog>
    lParalog = [ensembl_protein_ID for ensembl_protein_ID in my_gene_tree['leafs'].keys() if ensembl_protein_ID.startswith("ENSP0")]
    lParalog_Pfam = [iParalog for iParalog in lParalog if 'pfam repeats' in my_gene_tree['leafs'][iParalog]]

    # For each iPfam in <lPfam_ID>, find the min and max alignment coordinates of the corresponding TRs in all lParalog, that comply with lD>= 20 and nD >= 2.5
    # If at least one complying TR was found, create according <TR_file> in <results_dir>
    lTR_Pfam = [my_gene_tree['leafs'][iParalog]['pfam repeats'] for iParalog in lParalog if 'pfam repeats' in my_gene_tree['leafs'][iParalog]]
    lTR_Pfam = list(itertools.chain(*lTR_Pfam))

    for iPfam in lPfam_ID:
        lTR = [iTR for iTR in lTR_Pfam if iTR.Pfam == iPfam]
        if any(iTR.nD >= 2.5 and iTR.lD >= 10 and iTR.pValue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] < 0.1 for iTR in lTR):
            min_begin_alignment = min(iTR.position['alignment']['begin'] for iTR in lTR)
            max_end_alignment = max(iTR.position['alignment']['end'] for iTR in lTR)        
            homologs = {iParalog : [i for i in my_gene_tree['leafs'][iParalog]['pfam repeats'] if hasattr(i, 'Pfam') and i.Pfam == iPfam] for iParalog in lParalog_Pfam}
            homologs = {iParalog: jTR[0] for iParalog, jTR in homologs.items() if jTR != []}
            dHMM = {'PFAM_ID': iPfam}
            count_TR = create_tandem_repeat_dict(result_dir, my_gene_tree, count_TR, dHMM, homologs)
            coordinates_blocked.append((min_begin_alignment,max_end_alignment))
    
    # Find all TRs with l>=20, divergence <= 0.8, alpha <= 0.001 (gamma = 0.001), that are outside pfam boundaries.
    # Cluster TRs, and find best in cluster.
            
    # Find all tandem repeats in the current leaf that comply with the limitations and assign identifier
    tandem_repeats = []
    for iParalog in lParalog:
        iTandem_repeats = my_gene_tree['leafs'][iParalog]['tandem repeats']['aa'] 
        for i in iTandem_repeats:
            i.ensembl_ID = iParalog
            i.calc_nD()            
        tandem_repeats += [iTR for iTR in iTandem_repeats if not overlap(coordinates_blocked, [(iTR.position['aa_alignment']['begin'], iTR.position['aa_alignment']['end'])]) and iTR.lD >= c['lD'] and iTR.divergence[c['classifier']] <= c['d'] and iTR.pValue[c['classifier']] <= c['alpha']]
        
    if len(tandem_repeats) == 0:
        return count_TR
    
    # Cluster tandem repeats
    tandem_repeats = repeat_info.cluster_tandem_repeats(tandem_repeats, reference = 'aa_alignment')
    # Find best tandem repeat in each cluster
        
    for iCluster in tandem_repeats:
        # Find best tandem repeat in each cluster
        my_TR_seed = repeat_info.find_best_in_cluster(c,iCluster)
        
        # Build HMM from TR, and rerun on <lParalog>. If any single one complies with lD >= 20 and n >= 2.5, 
        # create according <TR_file> in <results_dir>
        my_hmm = hmm.HMM(tandem_repeat = my_TR_seed)
        
        homologs = {}
        for iParalog in lParalog: 
            leaf = my_gene_tree['leafs'][iParalog]
            protein_sequence = str(leaf['sequences']['aa'].seq).replace("*","").replace("X","")
            # Run Protein TRDs on Protein Sequence
            my_Viterbi = viterbi.Viterbi(my_hmm, protein_sequence)        
            my_Most_likely_path = my_Viterbi.viterbi()
            unaligned_msa = viterbi.hmm_path_to_non_aligned_tandem_repeat_units(protein_sequence, my_Most_likely_path, my_hmm.lD)
            if unaligned_msa == None or len(unaligned_msa) == 1:
                continue
            # Realign the repeat units
            aligned_msa = repeat_realign.realign_repeat(unaligned_msa)
            my_TR = repeat_info.Repeat(aligned_msa, calc_score = True, calc_pValue = True)
            my_TR.begin = my_Most_likely_path.count('N') + 1
        
            # Save start and end position in alignment.
            my_TR.position = {}
            my_TR.position['alignment'] = repeat_info.calculate_position_in_alignment(my_TR.begin, my_TR.sequence_length,alignment = leaf['sequences']['aa_gene_tree_alignment'])  
            my_TR.position['sequence'] = {'begin': my_TR.begin-1, 'end': my_TR.begin + my_TR.sequence_length - 1}
            homologs[iParalog] = my_TR
        if any((iTR.lD >= 10 and iTR.nD >=2.5 and iTR.pValue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] < 0.1) for iTR in homologs.values()):
            dHMM = {'tandem_repeat_seed': my_TR_seed}
            count_TR = create_tandem_repeat_dict(result_dir, my_gene_tree, count_TR, dHMM, homologs) 
    
    return count_TR
    
   
def overlap(a,b):
    
    for i,j in itertools.product(a,b):
        if (i[0] > j[0] and i[0] < j[1]) or (i[1] > j[0] and i[1] < j[1]) or (i[0] <= j[0] and i[1] >= j[1]):
            return True
    else:
        return False
    

def create_tandem_repeat_dict(result_dir, my_gene_tree, count, dHMM, homologs):
    
    ''' 
    save a tandem repeat dict
    construct a name for the result file from <my_gene_tree> and <count>
    save <dHMM> into the result file.
    
    Parameter example:
    dHMM = {'HMM': hmm_object, 'PFAM_ID' = 'PF00064', 'more_info': 17}
    
    '''
    
    tr_name = "{}_{}.pickle".format(my_gene_tree['gene_tree_ID'],count)
    result_file = os.path.join(result_dir,tr_name)
    
    tr = {'HMM': dHMM, 'gene_tree_ID': my_gene_tree['gene_tree_ID'], 'homologs': homologs}
    
    with open(result_file,'wb') as rh:
        pickle.dump(tr, rh)
    
    return count + 1
    
    
def extract_tandem_repeats_per_seq(gene_tree_file, gene_tree_index, result_dir, species = 'homo_sapiens', summary = None):
    
    ''' Find all none-overlapping TRs in the <my_gene_tree> dictionary, complying to <c> '''
    
    iCount = 0
    if summary == None:
        summary = {'l':[],'n':[],'divergence':[], 'ensembl_protein_ID': []}
    
    logging.debug("<gene_tree_file>: {0}".format(gene_tree_file))
    with open(gene_tree_file, 'rb') as result_filehandle:
        my_gene_tree = pickle.load(result_filehandle)
            
    # Define some kind of cut-off-constraints c for which repeats to use and which not:
    c = {'alpha' : 0.01, 'classifier' : 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions', 'l' : 20, 'd': 0.8} 
    # Define classifiers of interest:
    classifiers = ['phylo', 'phylo_gap01', 'phylo_gap001']
    
    species_found = False
    for ensembl_protein_ID,leaf in my_gene_tree['leafs'].items():
        if leaf['species'] == species:
            species_found = True
            # Find all tandem repeats in the current leaf that comply with the limitations and assign identifier
            tandem_repeats = [iTR for iTR in leaf['tandem repeats']['aa'] if iTR.lD >= c['l'] and iTR.divergence[c['classifier']] <= c['d'] and iTR.pValue[c['classifier']] <= c['alpha']]
            # Cluster tandem repeats
            tandem_repeats = repeat_info.cluster_tandem_repeats(tandem_repeats, reference = 'aa_protein_sequence')
            # Find best tandem repeat in each cluster

            for iCluster in tandem_repeats:
                count = str(len(summary['l']))
                iCount += 1
                # Find best tandem repeat in each cluster
                best = repeat_info.find_best_in_cluster(c,iCluster)
                result_file = os.path.join(result_dir, count + '.pickle')
                result_file_faa = os.path.join(result_dir, count + '.faa')
                # Write MSA of best none-overlapping tandem_repeat to fasta-file
                with open(result_file_faa, 'w') as f:
                    f.write(">{0} gene_tree:{1} begin:{2} l:{3} n:{4} pValue:{5} divergence:{6}\n".format(ensembl_protein_ID, str(gene_tree_index), str(best.begin), str(best.lD), str(best.n), str(best.pValue[c['classifier']]), str(best.divergence[c['classifier']])))
                    f.write("\n".join(best.msa))
                homolog_seed = {'gene_tree': gene_tree_index, 'original_seed': {'sequence': leaf['sequences']['aa'], 'tandem_repeat': best, 'Ensembl_Protein_ID': ensembl_protein_ID}}
                # Store all none-overlapping tandem_repeats
                with open(result_file, 'wb') as f:
                    pickle.dump(homolog_seed,f)
                summary['l'].append(best.lD)
                summary['n'].append(best.n)
                summary['divergence'].append(best.divergence[c['classifier']])
                summary['ensembl_protein_ID'].append(ensembl_protein_ID)

    return iCount,summary
    
## Write summary to tsv file:
# r = '/cluster/scratch_xl/public/eschaper/tandem_repeat_evolution/results/gene_tree_homologs_human_20_d08/gene_tree_homologs_human_20_d08_summary.tsv'
# with open(r,'w') as b:
#     b.write("\t".join(['ensembl_protein_ID', 'divergence', 'l', 'n']) + '\n')
#     for i,iE in enumerate(g['ensembl_protein_ID']):
#         b.write("\t".join([iE, str(g['divergence'][i]), str(g['l'][i]), str(g['n'][i])]) + '\n')    
    
     

def apply_hmm_per_tr(my_TR,result_file):

    ''' Try out several TRs with HMM of different parameter setting '''  
  
    my_TR_seed = my_TR['original_seed']['tandem_repeat']
    my_sequence = sequence.Sequence()
    my_sequence.sequence = str(my_TR['original_seed']['sequence'].seq).replace("*","").replace("X","")
  
    priors_divergence = [{'type': 'alpha', 'value': i} for i in [0.01,0.05,0.1]]    
    priors_indel_insertion = [{'type': 'adaptive', 'factor': i} for i in [1000,1000,1000]] # {'type': 'fixed', 'mu': 0.001, 'sigma_squared': 0.0002}
    #priors_divergence = [{'type': 'alpha', 'value': i} for i in [0.01]]    
    #priors_indel_insertion = [{'type': 'adaptive', 'factor': i} for i in [1000]] # {'type': 'fixed', 'mu': 0.001, 'sigma_squared': 0.0002}
 
    
    my_TR['seed_remake'] = []
    with open(result_file,'w') as f:
        scores = ' '.join(['pValue:{0} divergence:{1}'.format(my_TR_seed.pValue[i], my_TR_seed.divergence[i]) for i in SCORES])
        f.write(">{0} gene_tree:{1} begin:{2} l:{3} n:{4} {5}\n".format(my_TR['original_seed']['Ensembl_Protein_ID'], str(my_TR['gene_tree']), str(my_TR_seed.begin), str(my_TR_seed.lD), str(my_TR_seed.n), scores))
        f.write("\n".join(my_TR_seed.msa))
    
        # Define HMM parameter range
        for iPD, iPI in zip(priors_divergence, priors_indel_insertion):
            # Create and run HMM
            iHMM = hmm.HMM(my_TR_seed, prior_divergence=iPD, prior_indel_insertion=iPI)
            iViterbi = viterbi.Viterbi(iHMM, my_sequence.sequence)
            iMost_likely_path = iViterbi.viterbi()
            unaligned_msa = viterbi.hmm_path_to_tandem_repeat_units([my_sequence.sequence], [iMost_likely_path], my_TR_seed.lD)
            if unaligned_msa == None:
                continue
                my_TR['seed_remake'].append(None)
            elif len(unaligned_msa[0]) == 1:
                # Only one repeat unit was detected.
                f.write("\n\n>iPD:{0} iPI:{1} begin:{2} l:{3} n:{4}\n".format(str(iPD), str(iPI), str(iMost_likely_path.count('N')), len(unaligned_msa[0][0]), '1'))
                f.write(unaligned_msa[0][0])
                my_TR['seed_remake'].append(None)
                continue    
            # Realign the repeat units
            aligned_msa = repeat_realign.realign_repeat(unaligned_msa[0])
            iTR = repeat_info.Repeat(aligned_msa,scoreslist=SCORESLIST, calc_score = True, calc_pValue = True)
            iTR.begin = iMost_likely_path.count('N')
            # Calculate pValue & divergence
            scores = ' '.join(['pValue:{0} divergence:{1}'.format(iTR.pValue[i], iTR.divergence[i]) for i in SCORES])
            f.write("\n\n>iPD:{0} iPI:{1} begin:{2} l:{3} n:{4} {5}\n".format(str(iPD), str(iPI), str(iTR.begin), str(iTR.lD), str(iTR.n), scores))
            f.write("\n".join(iTR.msa))
            my_TR['seed_remake'].append(iTR)
    
    with open(tr_file, 'wb') as tr_filehandle:
        pickle.dump(my_TR,tr_filehandle)
  
def apply_hmm_per_tr_filter(my_TR, min_nD = 2, result_file = None):
    
    ''' Choose best TR, add nD
        If a result_file is given, save results
        else, return best hit. '''
        
    try:    
        lTR = [my_TR['original_seed']['tandem_repeat']] + my_TR['seed_remake']
    except:
        lTR = [my_TR['original_seed']['tandem_repeat']]
    for iTR in lTR:
        if iTR != None:
            iTR.calc_nD()
    lTR = [(i,iTR) for i,iTR in enumerate(lTR) if (iTR != None and iTR.nD >= min_nD and iTR.pValue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] < 0.001 and iTR.divergence['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] <= 0.8)]
    divergences = [iTR.divergence['phylo_gap001_ignore_coherent_deletions'] for i,iTR in lTR]
    if len(lTR) > 0:
        min_index = divergences.index(min(divergences))
        former_index = lTR[min_index][0]
        best_TR = lTR[min_index][1]
        if not result_file:
            return best_TR
        if former_index == 0:
            text = 'original seed;'
        elif former_index == 1:
            text = 'HMM;alpha:0.01'
        elif former_index == 2:
            text = 'HMM;alpha:0.05'
        elif former_index == 3:
            text = 'HMM;alpha:0.1'        
        with open(result_file,'w') as f:
            scores = ' '.join(['{0} pValue:{1} divergence:{2}'.format(i, best_TR.pValue[i], best_TR.divergence[i]) for i in ['phylo_gap001_ignore_coherent_deletions']])
            f.write(">{0} gene_tree:{1} {2} begin:{3} l:{4} nD:{5} {6}\n".format(my_TR['original_seed']['Ensembl_Protein_ID'], str(my_TR['gene_tree']), text, str(best_TR.begin), str(best_TR.lD), str(best_TR.nD)[:4], scores))
            f.write("\n".join(best_TR.msa))
      
        return former_index
    else:
        return None 


def load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir = None):
    
    tandem_repeat_file = os.path.join(tandem_repeat_dir, str(tandem_repeat_ID) + '.pickle')
    logging.debug("<tr_file>: {0}".format(tandem_repeat_file))
    with open(tandem_repeat_file, 'rb') as tr_filehandle:
        my_TR_data = pickle.load(tr_filehandle)
    
    if not gene_tree_dir:
        return my_TR_data
    else:
        gene_tree_file = os.path.join(gene_tree_dir, str(my_TR_data['gene_tree_ID']) + '.pickle') # WARNING! Formely "gene_tree"
        logging.debug("<gene_tree_file>: {0}".format(gene_tree_file))
        with open(gene_tree_file, 'rb') as gene_tree_filehandle:
            my_gene_tree = pickle.load(gene_tree_filehandle)
        return my_TR_data, my_gene_tree


def homologous_tandem_repeats_with_regard_to_pfam_HMM(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, result_dir):

    '''
        Load <my_TR_data> and <my_gene_tree> corresponding to <tandem_repeat_file> and <gene_tree_dir>
        Create <my_TR_data> HMM.
        Run HMM on all homologs in <my_gene_tree>.
        Save resulting MSAs as <my_TR_data['homologs']>.
    
    '''
    my_TR_data, my_gene_tree = load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)
        
    if 'PFAM_ID' in my_TR_data['HMM']:
        # Build Pfam HMM
        ### Define hmmer file
        hmm_file = os.path.join(DATAROOT2, "Pfam-A.hmm")
        my_HMM = hmm.HMM(hmm_file = hmm_file, accession=my_TR_data['HMM']['PFAM_ID'])
    else:
        my_TR_seed = my_TR_data['HMM']['tandem_repeat_seed']
        my_HMM = hmm.HMM(tandem_repeat = my_TR_seed)
    
    my_TR_data['lD'] = my_HMM.lD
    my_TR_data['homologs'] = predict_homologous_tandem_repeats_gene_tree(my_HMM, my_HMM.lD, my_gene_tree)
    result_file = os.path.join(result_dir,tandem_repeat_ID+'.pickle')
    with open(result_file,'wb') as rh:
        pickle.dump(my_TR_data,rh)

def homologous_tandem_repeats(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, cut_off_constraint=None):

    ''' 
    Get a <seed_TR> from <tr_file> - Either the original TR, or reruns of HMM on the same sequence.
    Predict Homologous tandem repeats to this <seed_TR> on the <gene_tree> (defined <tr_file>) in <gene_tree_dir>.
    Save the Viterbi-Paths of all predicted homologous tandem repeats.
    Save the MSAs of all predicted homologous tandem repeats as a fasta file.
    Realign all sequences with the won information on homologous tandem repeats using PrographMSA, and save the resulting MSA.
    
    '''
    my_TR_data, my_gene_tree = load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)

    # Extract TR seed
    my_seed_TR = apply_hmm_per_tr_filter(my_TR_data, min_nD = 2)
    if not my_seed_TR:
        # There is no seed TR in this probe.
        return None
    
    my_TR_data['final_seed'] = my_seed_TR  
    # Save prediction results
    with open(tr_file, 'wb') as tr_filehandle:
        pickle.dump(my_TR_data, tr_filehandle)
#     
#     
#     # Predict homologs
#     # Build a HMM out of the tandem repeat. 
#     prior_divergence = {'type': 'alpha', 'value': 0.1}
#     prior_indel_insertion = {'type': 'adaptive', 'factor': 1000}
#     my_HMM = hmm.HMM(my_TR, prior_divergence=prior_divergence, prior_indel_insertion=prior_indel_insertion)
# 
#     my_TR_data['homologs'] = predict_homologous_tandem_repeats_gene_tree(my_HMM,my_seed_TR.lD, my_gene_tree)
# 

def make_prograph_MSA(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir, result_file_alignment, result_file_fasta = None):

    my_TR_data, my_gene_tree = load_pickles(tandem_repeat_ID, tandem_repeat_dir, gene_tree_dir)

    # Check that the TR is indeed part of the alignment, make sure the index is correct.
    for ensembl_ID, data in my_TR_data['homologs'].items():
        sequence = str(my_gene_tree['leafs'][ensembl_ID]['sequences']['aa'].seq).replace("*","").replace("X","")
        if sequence[0] != 'M':
            sequence = 'M' + sequence 
        if 'MAFFT' in data['TR']:
            test = data['TR']['MAFFT'].repeat_in_sequence(sequence, save_original_msa = True)
            if not test:
                raise Exception("The tandem repeat is not part of the sequence.")
        else:
            data['TR']['HMM']['position']['begin'] =  sequence.find(data['TR']['HMM']['repeat_unit']) + 1   

    # Realign all sequences in gene tree with tandem repeat annotations in <result_file_treks>
    sequences = {ensembl_ID: str(leaf['sequences']['aa'].seq).replace("*","").replace("X","") for ensembl_ID, leaf in my_gene_tree['leafs'].items()}
    sequences = {ensembl_ID: iSeq if iSeq[0]=='M' else 'M'+iSeq for ensembl_ID, iSeq in sequences.items()}
    
    # Without Prograph <no_force_align> flag: use indices starting on 0 (but never a TR that actually starts on the first letter, if it is an 'M' :) ). Else: Use indices starting on 1.
    tandem_repeats = {ensembl_ID: ([data['TR']['MAFFT'].msa_original, data['TR']['MAFFT'].begin, data['TR']['MAFFT'].sequence_length] if 'MAFFT' in data['TR'] else [[data['TR']['HMM']['repeat_unit']], data['TR']['HMM']['position']['begin'], len(data['TR']['HMM']['repeat_unit'])]) for ensembl_ID, data in my_TR_data['homologs'].items()}
    prograph_msa = align.align(sequences, aligner = 'ProGraphMSA', sequence_type = 'AA', tandem_repeat_annotations = tandem_repeats)     
    
    my_TR_data['PrographMSA'] = {'MSA': prograph_msa}
    
    # Calculate position in alignment. tandem_repeats[id][1] = begin, tandem_repeats[id][2] = sequence_length
    position = {}
    for id,alignment in prograph_msa.items():
        position[id] = repeat_info.calculate_position_in_alignment(tandem_repeats[id][1], tandem_repeats[id][2], alignment)
        
    my_TR_data['PrographMSA']['position'] = position
        
    median_begin = int(np.floor(np.median([i['begin'] for i in position.values()])))
    median_end = int(np.ceil(np.median([i['end'] for i in position.values()])))
    
    my_TR_data['PrographMSA']['median_begin'] = median_begin
    my_TR_data['PrographMSA']['median_end'] = median_end
    
    # Save Prograph MSA
    align.save_msa_fasta(prograph_msa, result_file_alignment, tandem_repeat_position = {'begin': median_begin, 'end': median_end})
    
    # Save prediction results
    tr_file = os.path.join(tandem_repeat_dir,tandem_repeat_ID+'pickle')
    with open(tr_file, 'wb') as tr_filehandle:
        pickle.dump(my_TR_data, tr_filehandle)
    
#     # Print prediction results
#     with open(result_file_fasta, 'w') as f:
#         for ensembl_ID, data in my_TR_data['homologs'].items():
#             if 'MAFFT' in data['TR']:
#                 iTR = data['TR']['MAFFT']
#                 try: 
#                     scores = ' '.join(['{0} pValue:{1} divergence:{2}'.format(i, iTR.pValue[i], iTR.divergence[i]) for i in ['phylo_gap001_ignore_coherent_deletions']])
#                 except:
#                     scores = 'None'
#                 f.write(">{0} gene_tree:{1} begin:{2} l:{3} nD:{4} {5}\n".format(ensembl_ID, str(my_TR_data['gene_tree']), str(iTR.begin), str(iTR.lD), str(iTR.nD)[:4], scores))
#                 f.write("\n".join(iTR.msa)+"\n\n")
#             elif 'HMM' in data['TR']:
#                 # Only one repeat unit was detected.
#                 iTR = data['TR']['HMM']
#                 f.write(">>{0} gene_tree:{1} begin:{2} l:{3} nD:1\n".format(ensembl_ID, str(my_TR_data['gene_tree']), str(iTR['position']['begin']), str(len(iTR['repeat_unit']))))
#                 f.write(iTR['repeat_unit'] + '\n\n') 
    
    return 'Done'


def predict_homologous_tandem_repeats_gene_tree(my_HMM, lD, my_gene_tree):

    ''' Predict tandem repeats homologous to <my_TR> on all amino acid sequences within
        <my_gene_tree>.
        Save the resulting <viterbi_path>.
        Calculate the TR MSA as predicted by the HMM, and realign using MAFFT.
        Calculate for both TRs their position within an <alignment>.
        Calculate for the realigned MSA p-Values.
        
        Also return viterbi hits that do not span more than one repeat unit.
        
        Return a dict:
        {Ensembl_ID: {'viterbi_path': , 'TR': {'MAFFT': , 'HMM': } }}
        '''

    result = {}
    
    # Iterate through all leafs in the gene tree again.
    for ensembl_protein_ID,leaf in my_gene_tree['leafs'].items():
    
        result[ensembl_protein_ID] = {'TR':{}}
        # Copy the amino acid sequence for each leaf into the sequence.Sequence object <my_sequence>.
        my_sequence = str(leaf['sequences']['aa'].seq).replace("*","").replace("X","")
        
        # Create a viterbi.Viterbi object from the search sequence <my_sequence> and the HMM <my_HMM>.
        my_viterbi = viterbi.Viterbi(my_HMM, my_sequence)
        
        # Run the Viterbi algorithm.
        viterbi_path = my_viterbi.viterbi()
        begin = viterbi_path.count('N')
        
        # Convert the Viterbi path into repeat units.
        unaligned_msa = viterbi.hmm_path_to_non_aligned_tandem_repeat_units(my_sequence, viterbi_path, lD)
        
        if not unaligned_msa:
            continue
        
        elif len(unaligned_msa) == 1:
            # If there is only one repeat unit, save it as a string directly.
            position = {'alignment': repeat_info.calculate_position_in_alignment(begin + 1, len(unaligned_msa[0]),alignment = leaf['sequences']['aa_gene_tree_alignment']),
                        'sequence': {'begin': begin, 'end': begin + len(unaligned_msa[0]) - 1},
                         'begin': begin + 1}
            result[ensembl_protein_ID] = {'repeat_unit': unaligned_msa[0], 'position': position}
            
        else:
            # Realign the repeat units
            aligned_msa = repeat_realign.realign_repeat(unaligned_msa)
            if aligned_msa == []:
                print("Problem with MAFFT alignment.")
                print(unaligned_msa)
                print(aligned_msa)
            # Build a TR from the realignments
            iTR = repeat_info.Repeat(aligned_msa, begin=begin, calc_score=False, calc_pValue=False)
            # Save start and end position in alignment.
            iTR.position = {}
            iTR.position['alignment'] = repeat_info.calculate_position_in_alignment(iTR.begin, iTR.sequence_length,alignment = leaf['sequences']['aa_gene_tree_alignment'])  
            iTR.position['sequence'] = {'begin': iTR.begin-1, 'end': iTR.begin + iTR.sequence_length - 2}
            
            iTR.begin = begin
            # For completeness, save the viterbi_path as proposed by the HMM also.
            iTR.viterbi_path = viterbi_path
            
            result[ensembl_protein_ID] = iTR
            
    return result
 

if __name__ == "__main__":       
    # Read Gene Trees from Ensembl Compara emf nhx file (Get: Ensembl Protein&Gene ID, species, Tree Structure)
    # Fore each gene_tree from Ensembl Compara...
    #gene_tree_file = 'Compara.69.protein.nhx.emf_human_ribonuclease_inhibitor_lrr'
    gene_tree_id = 5
    gene_tree_file = 'Compara.69.protein.nhx.'+str(gene_tree_id)+'.emf'
    
    if 1 == 2:
        cut_off_constraint = {'classifier': 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions', 'alpha': 0.01}
        result_file = os.path.join(RESULTROOT,"gene_tree","0.pickle")
        with open(result_file, 'rb') as result_filehandle:
            gene_tree = pickle.load(result_filehandle)
        tandem_repeats = repeat_info.cluster_tandem_repeats(gene_tree,cut_off_constraint)
        for iTR in tandem_repeats:
            print(iTR)
            get_tandem_repeat_homologs(iTR, gene_tree, cut_off_constraint)
            
    
    if 1 == 3:
        count_tandem_repeat = 0
        with open(os.path.join(DATAROOT2,'Compara',gene_tree_file), 'r') as gene_tree_filehandle:
            for i,gene_tree in enumerate(tree_IO.get_gene_tree_info(gene_tree_filehandle)):
                count_tandem_repeat += assemble_gene_tree_data(gene_tree, str(i), 0)
                # At the moment, handle only one tree per time
                break

        
