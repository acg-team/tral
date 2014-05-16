import re

################################# TREE HOME BREW #########################################


class Sequence:
    
    def __init__(self, species, protein_id, gene_id):
        self.species = species
        self.protein_id = protein_id
        self.gene_id = gene_id

class Node:
    """ A phylogenetic Node """
    
    def __init__(self):
        self.name = 'Node1'
        self.children = []
        
    def get_child(self):
        self.children.append(Node())
        self.children[-1].set_mother(self)
        return self.children[-1]
    
    def set_mother(self, node):
        self.mother = node

class Tree:
    """ A phylogenetic tree """
    
    def __init__(self):
        
        self.name = 'Tree1'
        self.main_node = Node()

        

def read_ensembl_emf_nhx_line(line, sequence_data = None):
    
    ''' Parse an Ensembl NHX 
    Exemple execution: tree = read_ensembl_emf_nhx_line(line)
    Example read children info: tree.main_node.children[0].name            
    '''
    
    tree = Tree()
    current_mother = tree.main_node
    
    pat_new_child_node= re.compile("\(")
    pat_child_main_info = re.compile("([\w.]+):([\d.]+)")
    pat_child_additional_info = re.compile("\[&&NHX:(D=(?P<D>[YN]):)?(B=(?P<B>\d+):)?(G=(?P<G>[\w.]+):)?(T=(?P<T>\d+))?\]")
    
    pat_next_child = re.compile(",")
    pat_complete_child_node = re.compile("\)")
    
    pat_end = re.compile(";")
    
    count = 0
    old_length = 0
    
    while(line):
        if old_length == len(line):
            raise ValueError(line)   
        else:
            old_length = len(line)
        
        match = pat_new_child_node.match(line)
        if match:
            logger.log(LOGLEVEL_VERBOSE, match.group())
            current_mother = current_mother.get_child()
            count += 1
            line = line[len(match.group()):]
        
        match = pat_child_main_info.match(line)   
        if match:
            logger.log(LOGLEVEL_VERBOSE, match.group())
            line = line[len(match.group()):]
            current_mother.protein_ID = match.group(1)
            if sequence_data and current_mother.protein_ID in sequence_data:
                current_mother.species = sequence_data[current_mother.protein_ID]['species']
                current_mother.gene_ID = sequence_data[current_mother.protein_ID]['gene_ID']
            current_mother.distance = match.group(2)
        
        match = pat_child_additional_info.match(line)   
        if match:
            logger.log(LOGLEVEL_VERBOSE, match.group())
            line = line[len(match.group()):]
            
            # read additional info
            if match.groupdict()['B'] != None:
                current_mother.parent_branch_confidence_value = float(match.groupdict()['B'])
            if match.groupdict()['D'] != None:
                if match.groupdict()['D'] == 'T':
                    current_mother.duplication_node = True
                else:
                    current_mother.duplication_node = False
            if match.groupdict()['G'] != None:
                current_mother.gene_ID = match.groupdict()['G']
            if match.groupdict()['T'] != None:
                current_mother.taxonomy_ID = int(match.groupdict()['T'])
        
        match = pat_next_child.match(line)   
        if match:
            logger.log(LOGLEVEL_VERBOSE, match.group())
            line = line[len(match.group()):]
            # Add a sister
            current_mother = current_mother.mother.get_child()
            
        match = pat_complete_child_node.match(line)   
        if match:
            logger.log(LOGLEVEL_VERBOSE, match.group())
            line = line[len(match.group()):]
            # Return to mother
            current_mother = current_mother.mother
            
        match = pat_end.match(line)   
        if match:
            logger.log(LOGLEVEL_VERBOSE, "End OF LINE.")
            line = None
    
    return tree