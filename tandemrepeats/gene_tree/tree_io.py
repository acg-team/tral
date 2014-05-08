# Copyright (C) 2012 by Elke Schaper (elke@inf.ethz.ch) 

"""I/O function wrappers for the Ensembl emf nhx file format. 
   See: ftp://ftp.ensembl.org/pub/release-68/emf/ensembl-compara/ 
""" 

import logging, os, re, shutil, subprocess, tempfile
from . import align

from repeat.paths import *
from . import ensembl_species_info, repeat_tree_info


LOGLEVEL_VERBOSE = 8
logger = logging.getLogger(__name__)
logging.addLevelName(LOGLEVEL_VERBOSE, "VERBOSE") # for parser debug output

################################ ENSEMBL GENE TREE IO ####################################

def split_gene_tree_data(infile,outfile_path):

    ''' Split Ensembl Compara.69.protein.nhx.emf file in single gene tree files '''

    count = 0
    output_file = os.path.join(outfile_path,'Compara.69.protein.nhx.'+str(count)+'.emf')
    print(output_file)
    f = open(output_file,'w')
    for i, line in enumerate(infile):
        f.write(line)
        if line[:2] == "//":
            count += 1
            f.close()
            output_file = os.path.join(outfile_path,'Compara.69.protein.nhx.'+str(count)+'.emf')
            f = open(output_file,'w')
    f.close()
    print("I generated {0} gene_tree files".format(str(count)))
    

def get_gene_tree_info(infile, gene_tree_index = None): # Formerly known as read_ensembl_emf_nhx()
    """Read Trees from Ensembl emf nhx  files successively
    
    Arguments:
    infile: File stream
    
    Returns a generator function yielding one gene tree and all its sequences per request.
    
    If gene_tree_index is defined, only the gene_tree with the corresponding index is yielded.
    
    Postcondition: infile points to EOF.
    """
    
    # pattern for start index
    pat_seq_data = re.compile("SEQ (?P<species>[\w_]+) (?P<protein_ID>[\w.]+) [\S]+ [\d]+ [\d]+ [\d-]+ (?P<gene_ID>[\w.]+)")
    pat_data = re.compile(r"DATA")
    
    # Our possible parser states:
    #
    # 1: searching for sequence data
    # 2: saving tree
    
    if gene_tree_index:
        count = 0
    
    sequence_data = {}
    state = 1
    for i, line in enumerate(infile):
        logger.log(LOGLEVEL_VERBOSE, "Line %d: %s", i, line[0:-1])
        if 1 == state:
            match = pat_seq_data.search(line)
            if match:
                logger.log(LOGLEVEL_VERBOSE, " * (1->1) Found sequence data")
                sequence_data[match.groupdict()['protein_ID']] = {'species': match.groupdict()['species'], 'gene_ID': match.groupdict()['gene_ID']}
                continue
            elif pat_data.search(line):
                state = 2
                logger.log(LOGLEVEL_VERBOSE, " * (1->2) Found data tag")
                continue
        
        if 2 == state:
            state = 1
            logger.log(LOGLEVEL_VERBOSE, " * (2->1) Found a tree.")
            ## Execute this line, if you would like to yield a Tree Instance, and not a string
            #yield read_ensembl_emf_nhx_line(line, sequence_data)
            if gene_tree_index == None or count == gene_tree_index:
                yield {'leafs': sequence_data, 'tree_nhx': line}
            count = count + 1
            continue
            # save



################################## TREE VISUALISATION ####################################
    
def create_scriptree_files(result_file_stump, tree = None,  tree_nhx = None, annotations = None):

    if not tree:
        tree = repeat_tree_info.Repeat_Tree(tree_nhx)
    if not annotations:
        # Create Tree annotations
        annotations = {}  
        for iLeaf in tree.iter_leaves():
            iEnsembl_ID, iUnit = iLeaf.name.split("_")
            species = ensembl_species_info.translate('ensembl_ID','official', iEnsembl_ID)
            annotations[iLeaf.name] = 'Species {{{0}}} EnsemblID {{{1}}} RepeatUnit {{{2}}}'.format(species, iEnsembl_ID, iUnit)
    
    annotation_file = os.path.join(result_file_stump + '_annotations.txt')
    with open(annotation_file, 'w') as annotation_handle:
        for id, data in annotations.items():
            annotation_handle.write("{0} {1}\n".format(id, data))
            
    tree.write(format=5, outfile=os.path.join(result_file_stump + '.nhx'))

def draw_tree(tree_file,result_file, annotation_file = None, script_file = None, xvfb = False):

    '''
      Create ps files from trees using script tree: http://lamarck.lirmm.fr/scriptree/
      Convert the .ps to .pdfs?
    '''
    
    if not annotation_file:
        annotation_file = os.path.join(EXECROOT, 'scriptree','annotations.txt')
        
    if not script_file:
        script_file = os.path.join(EXECROOT, 'scriptree','script.txt')
    scripttree_dir = os.path.join(EXECROOT, 'scriptree')
    scriptree_bin = os.path.join(scripttree_dir,'scriptree.bin')
    #scriptree = os.path.join(scripttree_dir,'scriptree.tcl')
    #etcl = os.path.join(EXECROOT,'etcl/bin/etcl')
    
    commands = [scriptree_bin, "-tree",  tree_file, "-annotation",  annotation_file, "-script", script_file, "-out", result_file, "2>/dev/null"]
    if xvfb:
        commands = ['/cluster/apps/xvfb-run/xvfb-run'] + commands

    logging.debug(" ".join(commands))
    try:
        p = subprocess.Popen(commands, cwd = scripttree_dir)
        p.wait()
    except:
        return commands, scripttree_dir
    

if __name__ == "__main__":       

    gene_tree_file = 'Compara.69.protein.nhx.emf'
    from repeat.paths import *
    if 1 == 1:
        gene_tree_file = os.path.join(DATAROOT2,gene_tree_file)
        with open(gene_tree_file, 'r') as filehandle:
            split_gene_tree_data(filehandle,os.path.join(DATAROOT2,"Compara"))
