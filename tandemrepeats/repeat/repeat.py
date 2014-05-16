# (C) 2011, Alexander Korsunsky
# (C) 2011-2013 Elke Schaper

from os.path import join
import re
import logging
import numpy as np
import scipy.special
from collections import defaultdict
from copy import deepcopy

from . import repeat_score, repeat_pvalue

logger = logging.getLogger(__name__)

PRECISION = 10
SCORESLIST = ['phylo', 'phylo_gap01', 'phylo_gap001']

################################### Repeat class #########################################



class Repeat:
    
    def __str__(self):
        if hasattr(self, 'pValue'):
            first_line = '>begin:{0} lD:{1} n:{2} pValue:{3} divergence:{4}\n'.format(self.begin, self.lD, self.n, self.pValue[SCORESLIST[-1]], self.divergence[SCORESLIST[-1]])
        else:
            first_line = '>begin:{0} lD:{1} n:{2}\n'.format(self.begin, self.lD, self.n)
    
        return first_line + "\n".join(self.msa)

    def __init__(self, msa, begin = None, job_prefix="",
                 sequence_type = 'AA',
                 scoreslist=SCORESLIST,
                 calc_score = False, calc_pValue = False):
        
        # The index of begin start out on Zero.
        if begin != None:
            self.begin = begin
        
        # Replace underscores and asterisks with scores, transform letters to uppercase, replace Selenocysteine (U) with Cysteine (C), replace Pyrrolysine (O) with Lysine (K)
        self.msa_original = [repeat_unit.upper().replace("_", "-").replace("*", "-") for repeat_unit in msa]
        self.msa = [repeat_unit.replace("U", "C").replace("O", "K") for repeat_unit in self.msa_original]
        self.sequence_type = sequence_type

        self.n = len(self.msa)
        self.l = max(map(len, self.msa))

        # Assure every line is equally long
        for i, rep in enumerate(self.msa):
            self.msa[i] += "-" * (self.l-len(msa[i]))
        logging.debug("Constructing repeat from MSA:\n%s" % "\n".join(self.msa))

        ## Transpose MSA
        self.msaT = ["".join(c) for c in zip(*self.msa)]

        # Get text for MSA.
        self.text = " ".join(self.msa)

        self.nGap = self.text.count("-")
        self.sequence_length = len(self.text) - self.nGap - self.text.count(" ")

        self.deleteInsertionColumns()
        if self.lD != 0:
            self.gapStructure()

            if calc_score:
                self.calculate_scores(scoreslist=scoreslist, job_prefix=job_prefix)
    
            if calc_pValue:
                self.calculate_pValues(scoreslist=scoreslist)
         
    def calc_nD(self):
        ''' what is the number of letters in lD, divided by lD?'''  
        
        if self.lD != 0:
            self.nD = (len("".join(self.msaD)) - self.textD.count('-'))/self.lD
        else:
            self.nD = 0
        
    def calc_index_msa(self):
        ''' formerly named: set_begin'''

        ## calculate the transposed MSA with the index in protein sequence instead of letter
        ## to be able to compare two overlapping Repeats within the same sequence region.
        index = list(range(self.begin + self.sequence_length -1, self.begin -1, -1)) ## python 2: range returns list by defaul
        # Replace every letter with the index in the protein sequence
        msaI = [[j if j == '-' else index.pop() for j in i] for i in self.msa] ## the 'I' stands for index
        ## Transpose msaI using list-comprehensions/zip/join magic
        self.msaIT = [[i for i in c if i != '-'] for c in zip(*msaI)]     
        self.msaIT = [i for i in self.msaIT if len(i) > 0]

    def calculate_scores(self, scoreslist=SCORESLIST,
        job_prefix=""):

        if not hasattr(self,'score'):
            self.score = defaultdict(int)
            
        if not hasattr(self,'divergence'):
            self.divergence = defaultdict(int)
        
        worst_score = {
            'entropy' : 1,
            'parsimony' : 1,
            'pSim' : 0,
            'phylo' : 0,
            'phylo_gap01': 0,
            'phylo_gap001': 0}

        if  self.lD == 0: # Enter worst score
            self.score = {iScore: worst_score[iScore] for iScore in scoreslist}   
        elif self.n < 2: # Enter None, as a score cannot be calculated if there is just one repeat unit.
            self.score = {iScore: None for iScore in scoreslist}  
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.score['entropy'] = repeat_score.meanSimilarity(self, repeat_score.entropy)
            if 'parsimony' in scoreslist:
                self.score['parsimony'] = np.round_(repeat_score.meanSimilarity(self, repeat_score.parsimony), decimals = PRECISION)
            if 'pSim' in scoreslist:
                self.score['pSim'] = np.round_(repeat_score.meanSimilarity(self, repeat_score.pSim), decimals = PRECISION)
            if 'phylo' in scoreslist:
                self.divergence['phylo'], self.score['phylo'] = repeat_score.phyloStarTopology_local(self) 
            if 'phylo_gap01' in scoreslist:
                self.divergence['phylo_gap01'], self.score['phylo_gap01'] = repeat_score.phyloStarTopology_local(self, gaps = 'row_wise', indelRatePerSite = 0.01)
                self.divergence['phylo_gap01_ignore_trailing_gaps'], self.score['phylo_gap01_ignore_trailing_gaps'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps', indelRatePerSite = 0.01)
                self.divergence['phylo_gap01_ignore_coherent_deletions'], self.score['phylo_gap01_ignore_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_coherent_deletions', indelRatePerSite = 0.01)
                self.divergence['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'], self.score['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps_and_coherent_deletions', indelRatePerSite = 0.01)
            if 'phylo_gap001' in scoreslist:
#                self.divergence['phylo_gap001'], self.score['phylo_gap001'] = repeat_score.phyloStarTopology_local(self, gaps = 'row_wise', indelRatePerSite = 0.001)
#                self.divergence['phylo_gap001_ignore_trailing_gaps'], self.score['phylo_gap001_ignore_trailing_gaps'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps', indelRatePerSite = 0.001)
#                self.divergence['phylo_gap001_ignore_coherent_deletions'], self.score['phylo_gap001_ignore_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_coherent_deletions', indelRatePerSite = 0.001)
                self.divergence['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'], self.score['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps_and_coherent_deletions', indelRatePerSite = 0.001)
    
    def calculate_pValues(self, scoreslist=SCORESLIST):
       
        if not hasattr(self,'pValue'):
            self.pValue = defaultdict(int) 

        if  self.lD == 0: # Enter worst p-value
            self.pValue = {iScore: 1.0 for iScore in scoreslist} 
        elif self.n < 2: # Enter None, as a pValue cannot be calculated if there is just one repeat unit.
            self.pValue = {iScore: None for iScore in scoreslist}           
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.pValue['entropy'] = repeat_pvalue.pValueFromEmpiricialList(self, 'entropy')
            if 'parsimony' in scoreslist:
                self.pValue['parsimony'] = repeat_pvalue.pValuePars(self)
            if 'pSim' in scoreslist:
                self.pValue['pSim'] = repeat_pvalue.pValuePSim(self)
            if 'phylo' in scoreslist:
                self.pValue['phylo'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo')              
            if 'phylo_gap' in scoreslist:
                self.pValue['phylo_gap'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap')                
            if 'phylo_gap01' in scoreslist:
                self.pValue['phylo_gap01'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01')
                self.pValue['phylo_gap01_ignore_trailing_gaps'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score['phylo_gap01_ignore_trailing_gaps'])
                self.pValue['phylo_gap01_ignore_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score['phylo_gap01_ignore_coherent_deletions'])
                self.pValue['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'])
            if 'phylo_gap001' in scoreslist:
#                self.pValue['phylo_gap001'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001')
#                self.pValue['phylo_gap001_ignore_trailing_gaps'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score['phylo_gap001_ignore_trailing_gaps'])
#                self.pValue['phylo_gap001_ignore_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score['phylo_gap001_ignore_coherent_deletions'])
                self.pValue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'])
                   

    def deleteInsertionColumns(self):
        ''' Create the MSA and the transposed MSA for this tandem repeat lacking insertion columns.'''
                
        ## By definition, insertion columns are columns were there are more or
        ## equal gaps compared to amino acids.
        self.msaTD = [column
            for column in self.msaT if (2*column.count('-') < len(column))
        ]
        self.insertions = re.compile("-+").findall(''.join([
                'p' if (2*column.count('-') < len(column)) else '-'
                for column in self.msaT
            ])
        )

        self.msaD = ["".join(c) for c in zip(*self.msaTD)]

        self.lD = len(self.msaTD)
        self.textD = " ".join(self.msaD)
        self.totD = len(self.textD) - self.textD.count('-') - self.n + 1
        if self.lD != 0:
            self.nD = (len("".join(self.msaD)) - self.textD.count('-'))/self.lD
        else:
            self.nD = 0

    def gapStructure(self):
        ''' Calculate the number and length of insertions and deletions for this tandem repeat.'''
        
        ## Needed: Number and size of insertions, number and size of deletions
        # row_wise: (standard case) column wise coherence is ignored, the indels are summarised over the columns.
        # ignore_trailing_gaps: In addition to the standard case, trailing gaps on either side of the tandem repeat are replaced with letters, and hence ignored.
        # The separation into insertions and deletions comes before ignoring trailing gaps. 
        # ignore_coherent_deletions : In addition to the standard case, same size deletions spanning the same columns are only counted once.
        self.insertions = []
        self.deletions = {}
        self.gaps = {}
    
        ## row_wise
        # 1. Detect insertions
        # First, name columns either '-' if they are insertion columns, or 'p', if they are not.
        # Concetanate the string of 'p's and '-', and double it, to consider insertions in the last column
        # and the first column as associated (= of the same origin)
        # Lastly, detect the length of all insertions.
        insertions = re.compile("-+").findall(2*(''.join([
                'p' if (2*column.count('-') < len(column)) else '-' 
                for column in self.msaT
                ])))
        insertions = insertions[int(len(insertions)/2):-1] if len(insertions)%2 == 1 else insertions[int(len(insertions)/2):]
        self.insertions = [len(i) for i in insertions]
        if self.lD == 0:
            self.insertions = [self.l]
               
        # 2. Detect deletions
        # CHECK this lines. you used msaTD before, but swapped to msaD
        deletions = [(m.start()%self.l, len(m.group())) for m in re.finditer(re.compile("-+"), "".join(self.msaD))]
        self.deletions['row_wise'] = [i[1] for i in deletions]    
        self.gaps['row_wise'] = self.insertions + self.deletions['row_wise']
 
        ## ignore_coherent_deletions
        # Remove duplicates from the list of tuples(index of deletion, deletion length)
        deletions_ignore_coherent = {iD: None for iD in deletions} 
        self.deletions['ignore_coherent_deletions'] = [i[1] for i in deletions_ignore_coherent.keys()]
        self.gaps['ignore_coherent_deletions'] = self.insertions + self.deletions['ignore_coherent_deletions']
    
        ## ignore_trailing_gaps
        # To be precise, we should ignore_trailing_gaps only, if we have ascertained that there were no letters 
        # in possible insertion columns before the first deletion column. For now, this case is ignored.
        msaD = deepcopy(self.msaD)
        def sub_trailing_gap(match_object):
            return('x'*len(match_object.group())) 
        msaD[0] = re.sub('^-+', sub_trailing_gap, msaD[0])
        msaD[-1] = re.sub('-+$', sub_trailing_gap, msaD[-1])
        
        # CHECK this lines. you used msaTD before, but swapped to msaD
        #msaTD = ["".join(c) for c in zip(*msaD)]
        deletions = [(m.start()%self.l, len(m.group())) for m in re.finditer(re.compile("-+"), "".join(msaD))]
        self.deletions['ignore_trailing_gaps'] = [i[1] for i in deletions]
        self.gaps['ignore_trailing_gaps'] = self.insertions + self.deletions['ignore_trailing_gaps']
        
        ## ignore_trailing_gaps & ignore_coherent_deletions
        deletions_ignore_coherent = {iD: None for iD in deletions} 
        self.deletions['ignore_trailing_gaps_and_coherent_deletions'] = [i[1] for i in deletions_ignore_coherent.keys()]
        self.gaps['ignore_trailing_gaps_and_coherent_deletions'] = self.insertions + self.deletions['ignore_trailing_gaps_and_coherent_deletions']

    def gap_structure_HMM(self):

        ''' calculate for each column in self.msaTD 
        - how many deletions of what length begin in this column?: <deletion_lengths>[<iColumn>] = list of deletion lengths
        - how many insertions of what length begin in potential insertion columns before this column?: <insertion_lengths>[<iColumn>] = list of insertion lengths
        '''

        ## Insertions: In which column (with reference to self.msaTD) do how many insertions of what length start?
        # In which columns are insertions?
        insertions = [i for i,column in enumerate(self.msaT) if (2*column.count('-') >= len(column))]
        
        # Which insertions are neighboring == form a block, and on which index with reference to self.msaTD do they start?
        insertion_blocks = {}
        nInsertion_columns = 0  
        while insertions:
            current = insertions.pop(0)
            insertion_blocks[current] = {'indices_self.msaT':[current]}
            while insertions and insertions[0] == insertion_blocks[current]['indices_self.msaT'][-1] + 1:
                insertion_blocks[current]['indices_self.msaT'].append(insertions.pop(0))
            insertion_blocks[current]['index_self.msaTD'] = insertion_blocks[current]['indices_self.msaT'][0] - nInsertion_columns
            nInsertion_columns += len(insertion_blocks[current]['indices_self.msaT'])
        
        # If an insertions ranges over the repeat unit border, take the insertion on both parts together:
        if 0 in insertion_blocks.keys() and len(self.msaT)-1 in insertion_blocks.keys():
            insertion_blocks[current]['indices_self.msaT'].extend(insertion_blocks[0]['indices_self.msaT'])
            insertion_blocks[current]['index_self.msaTD'] = insertion_blocks[0]['index_self.msaTD']
            del insertion_blocks[0]
        
        # Restructure insertion_blocks: We are only interested in the index with reference to self.msaTD, and the start and end index with reference to self.msaT:
        insertion_blocks = {i['index_self.msaTD']: (i['indices_self.msaT'][0],i['indices_self.msaT'][-1]) for i in insertion_blocks.values()}
        
        # Go to each insertion block, and see how many insertion it includes, and of what length.
        l = len(self.msaT)
        insertion_lengths = {}
        for iL in range(self.lD):
            insertion_lengths[iL] = []
            if iL in insertion_blocks:
                begin = insertion_blocks[iL][0]
                end = insertion_blocks[iL][1]
                length = end - begin + 1
                # In case the insertion blocks bridges over then end of the repeat unit:
                if begin > end:
                    length += l
                print(length)
                    
                repeat_sequence = ''.join(self.msa)
                for i in range(self.n):
                    insertion = repeat_sequence[i*l + begin:min(len(repeat_sequence),i*l + begin+length)]
                    insertion_lengths[iL].extend([len(iI) for iI in re.findall('\w+',insertion)])
                
                # In case the insertion blocks bridges over then end of the repeat unit:
                if begin > end:
                    insertion = repeat_sequence[:end+1]
                    insertion_lengths[iL].extend([len(iI) for iI in re.findall('\w+',insertion)])
           
        logging.debug("The insertion_lengths of the current TR are: {0}".format(insertion_lengths))
                 
        ## Deletions are much easier to handle than insertions, as they are already in the right coordinate system (self.msaTD)
        # Find all deletions, and note their start index, and their length
        deletions_all = [(m.start(),len(m.group())) for m in re.finditer('-+', ''.join(self.msaD))]
        # Group the deletions with respect to the column in which they start:
        deletion_lengths = {iL: [iD[1] for iD in deletions_all if iD[0]%self.lD == iL] for iL in range(self.lD)}
        logging.debug("The deletion_lengths of the current TR are: {0}".format(deletion_lengths))
        
        return insertion_lengths, deletion_lengths

    def isSameRepeat(self, repeat2):
        ''' return 1 if the two TRs share at least one pair of amino acids (or
        nucleotides) with common ancestry;
        else 0 '''
        if not hasattr(self,'msaIT'):
            self.calc_index_msa()
        if not hasattr(repeat2, 'msaIT'):
            repeat2.calc_index_msa()   
        
        original = repeat2.msaIT
        potential = self.msaIT
        try:
            for p in potential:
                iP = 0
                for iP in range(len(p)-1):
                    i = 0
                    j = 0
                    while original[i][j] != p[iP]:
                        if original[i][j] > p[iP]:
                            i += 1
                            j = 0
                        else:
                            j += 1
                        if len(original) <= i or len(original[i]) <= j:
                            break
                    else:
                        for iPRest in range(iP + 1, len(p)):
                            if p[iPRest] in original[i][j+1:]:
                               self.isSame = True
                               self.coverage = self.sequence_length/repeat2.sequence_length
                               self.greediness = self.lD/repeat2.lD
                               return True
        except:
            logging.warning('error in isSameRepeat with original %s and potential %s',
                str(self.msaIT), str(repeat2.msaIT))
            logging.warning('original:', str(repeat2.msa))
            logging.warning('original:', str(repeat2.begin))
            logging.warning('potential:', str(self.msa))
            logging.warning('potential:', str(self.begin))
            self.isSame = False
            return False
        self.isSame = False
        return False
    
    def repeat_in_sequence(self,sequence, save_original_msa = False):
    
        ''' SANITY CHECK whether the predicted TR truly is part of the string <sequence>
        in which it was detected.
        
        If yes: Return True, set self.begin to corrected value if necessary.
        If no: Return False. '''

        repeat_sequence = "".join(self.msa).upper().replace("_", "").replace("-", "")
        starts = [m.start()+1 for m in re.finditer(repeat_sequence,sequence.replace("U", "C").replace("O", "K"))] # The first letter in the sequence is counted as 1 (not 0, as in python).
                
        if self.begin in starts: # Was the begin predicted correctly?
            self.save_original_msa(sequence)
            return True
        elif len(starts) != 0: # is the tandem repeat somewhere else?
            self.begin = starts[0]
            self.save_original_msa(sequence)
            return True
        else:
            return False


    def save_original_msa(self, sequence):
    
        ''' recalculate the original msa, assuming that begin is defined correctly (index starting counting on 1.) '''
        repeat_sequence = sequence[self.begin-1:self.begin-1 + self.sequence_length]
        
        count = 0
        self.msa_original = []
        for unit in self.msa:
            unit_original = ''
            for char in unit:
                if char == '-':
                    unit_original += '-'
                else:
                    unit_original += repeat_sequence[count]
                    count += 1
            self.msa_original.append(unit_original)
            
################################## Repeat functions ######################################
        
            
def calculate_position_in_alignment(begin, length, alignment):

    ''' calculate the index of the begin and the end of a TR within an alignment
        returns  
     '''
    
    # a alignment.seq (hopefully!) is a string of letters, with many gaps "-", as it this particular case is an alignment sequence.
    seq_indexed = [i for i,j in enumerate(str(alignment)) if j != '-']
    logger.debug("begin: {0}; length: {1}".format(str(begin), str(length)))
    logger.debug("alignment: {0}".format(str(alignment)))
    return {'begin': seq_indexed[begin-1], 'end': seq_indexed[begin + length - 2]}
                

def cluster_tandem_repeats(tandem_repeats, reference = 'aa_protein_sequence'):
    ''' (At current) Return best tandem repeat in each tandem repeat cluster.
    A cluster contains all overlapping tandem repeat with <reference> to 'aa_alignment'
    or 'aa_protein_sequence'. '''

    # Cluster these according to which part of the alignment they cover. Anything that overlaps ends up in the same cluster.    
    clusters = []
    count = 0
    while(tandem_repeats):
        cluster = [tandem_repeats.pop()]
        cluster_begin = cluster[0].position[reference]['begin'] 
        cluster_end = cluster[0].position[reference]['end'] 
        complete_cluster, tandem_repeats = find_overlapping_tandem_repeats(cluster_begin,cluster_end,tandem_repeats, cluster, reference)
        for iTR in complete_cluster:
            iTR.cluster = {'position_cluster': count}
        count += 1
        clusters.append(complete_cluster)

    logging.debug("{0} clusters were found.".format(str(len(clusters))))
    return clusters
    
def find_overlapping_tandem_repeats(cluster_begin,cluster_end,tandem_repeats, cluster, reference):
    ''' find all tandem_repeats overlapping cluster_begin or cluster_end. 
        (recursive, as cluster_begin and cluster_end can be changed dynamically) '''

    remaining_tandem_repeats = tandem_repeats[:]
    for i,iTR in enumerate(tandem_repeats):
        # Is tandem repeat within the boundaries?
        if iTR.position[reference]['begin'] >= cluster_begin and iTR.position[reference]['end'] <= cluster_end:
            cluster.append(iTR)
            remaining_tandem_repeats.remove(iTR)
            #cluster.append(remaining_tandem_repeats.pop(i)) 
        # Is the tandem repeat outside the boundaries?
        elif iTR.position[reference]['end'] < cluster_begin or iTR.position[reference]['begin'] > cluster_end:
            continue
        else:    
            # Is the tandem_repeat spanning the current boundaries?    
            if iTR.position[reference]['begin'] < cluster_begin and iTR.position[reference]['end'] > cluster_end: 
                cluster_begin = iTR.position[reference]['begin']   
                cluster_end = iTR.position[reference]['end']
            # Is the tandem repeat starting a bit before the boundaries, but ending within?
            elif iTR.position[reference]['begin'] < cluster_begin and iTR.position[reference]['end'] <= cluster_end:
                cluster_begin = iTR.position[reference]['begin']
            # Is the tandem repeat ending a bit after the boundaries, but starting within?
            elif iTR.position[reference]['begin'] >= cluster_begin and iTR.position[reference]['end'] > cluster_end:                
                cluster_end = iTR.position[reference]['end']
            else:
                print("WARNING! IS THIS CASE ALLOWED TO HAPPEN?")
            
            # In all three cases, add the tandem_repeat to the cluster, and start a new search with new boundaries recursively.    
            cluster.append(iTR)
            remaining_tandem_repeats.remove(iTR)
            #cluster.append(remaining_tandem_repeats.pop(i))
            return find_overlapping_tandem_repeats(cluster_begin,cluster_end,remaining_tandem_repeats, cluster, reference)
    
    return cluster, remaining_tandem_repeats


def find_best_in_cluster(cut_off_constraint,cluster):

    '''
    Parameters:
    e.g.
    cut_off_constraint = {'classifier': 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions', 'alpha': 0.01}
    '''

    pValues = [iTR.pValue[cut_off_constraint['classifier']] for iTR in cluster]
    min_pValue = min(pValues)
    cluster = [cluster[i] for i,j in enumerate(pValues) if j == min_pValue]
    if len(cluster) > 1:
        logging.debug("{0} best tandem repeats had the same pValue.".format(str(len(cluster))))
        divergences = [iTR.divergence[cut_off_constraint['classifier']] for iTR in cluster]
        min_divergence = min(divergences)
        cluster = [cluster[i] for i,j in enumerate(divergences) if j == min_divergence]
        if len(cluster) > 1:
            logging.debug("{0} best tandem repeats had the same divergence.".format(str(len(cluster))))
            sequence_lengths = [iTR.sequence_length for iTR in cluster]
            max_sequence_length = max(sequence_lengths)
            cluster = [cluster[i] for i,j in enumerate(sequence_lengths) if j == max_sequence_length]
            if len(cluster) > 1:
                logging.debug("{0} best tandem repeats had the same maximum sequence length.".format(str(len(cluster))))
                lDs = [iTR.lD for iTR in cluster]
                min_lD = min(lDs)
                cluster = [cluster[i] for i,j in enumerate(lDs) if j == min_lD]
                if len(cluster) > 1:
                    logging.debug("{0} best tandem repeats had the same lD.".format(str(len(cluster))))
                    logging.debug("{0} best tandem repeats could not be destinguished. One is picked randomly.".format(str(len(cluster))))  
    return cluster[0]


    
def mark_best_in_clusters(cut_off_constraint,clusters):    

    # For each cluster: 
    # Witch TR has the smallest p-Value? If undecided, which has the smallest divergence? If undecided, which has the smallest lD? If undecided, which one is longest? If undecided: First in list.
    for iCluster in clusters:
        best_TR = find_best_in_cluster(cut_off_constraint,iCluster)
        # mark the best tandem repeat        
        best_TR.cluster['best'] = True 


########################################### MAIN #########################################   
        
if __name__=="__main__":

    myTR = Repeat(begin = 10, msa = ['DFG','DFG','DFG'], sequence_type = 'AA', calc_score = True, calc_pValue = True)
    

