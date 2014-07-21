# (C) 2011, Alexander Korsunsky
# (C) 2011-2014 Elke Schaper

from collections import defaultdict
import configobj
from copy import deepcopy
import logging
import numpy as np
import os
from os.path import join
import re
import scipy.special

from tandemrepeats.repeat import repeat_score, repeat_pvalue, repeat_io
from tandemrepeats.paths import *


logger = logging.getLogger(__name__)

pDefaults = os.path.join(CODEROOT, 'tandemrepeats', 'data', 'defaults.ini')
pSpec = os.path.join(CODEROOT, 'tandemrepeats', 'data', 'spec.ini')
configs = configobj.ConfigObj(pDefaults, configspec = pSpec)
config = configs["repeat"]

################################### Repeat class #########################################

class Repeat:

    """ Tandem repeats are consecutive copies of genomic sequence.

    A tandem repeat is stored in terms of a multiple sequence alignment of its
    tandem repeat units. For example, the amino acid tandem repeat
    KLMKLMKLKLM
    can be displayed as a multiple sequence alignment (here: "multiple repeat unit alignment")
    KLM
    KLM
    KL-
    KLM
    and can be stored as
    ["KLM", "KLM", "KL-", "KLM"].
    The - char stands for a gap: It indicates that either, there used to be a char at the
    same position, but it got deleted. Or, all other chars in the same column got inserted
    at some point.

    For research on tandem repeats, many different statistics need to be calculated on
    the multiple sequence alignment of the tandem repeat.

    [Description of these statistics]

    Attributes:
        msa_original (list of str): The original alignment of repeat units.
        msa (list of str):  The alignment of repeat units without rare chars.
        msaT (list of str): The transposed alignment of repeat units. E.g.
            ["KKKK", "LLLL", "MM-M"].
        sequence_type (str): "AA" or "DNA"
        n (int): The number of repeat units in the original msa. E.g. 4.
        l (int): The length of the repeat units in the original msa. E.g. 3.
        nD (int):  [Needs explanation]
        lD (float): [Needs explanation]
        text (str): [Needs explanation]
        nGap (int): The number of gaps in ``msa``
        sequence_length (int): [Needs explanation]
    """


    def __str__(self):

        """
        .. todo:: Improve output if pValues are only partly calculated.

        """

        if hasattr(self, 'pValue') and hasattr(self, 'begin'):
            first_line = '>begin:{0} lD:{1} n:{2} pValue:{3} divergence:{4}\n'.format(self.begin, self.lD, self.n, self.dPValue[config['scoreslist'][-1]], self.dDivergence[config['scoreslist'][-1]])
        elif hasattr(self, 'pValue'):
            first_line = '>lD:{} n:{} pValue:{} divergence:{}\n'.format(self.lD, self.n, self.dPValue[config['scoreslist'][-1]], self.dDivergence[config['scoreslist'][-1]])
        else:
            first_line = '>lD:{} n:{}\n'.format(self.lD, self.n)

        return first_line + "\n".join(self.msa)

    def write(self, format, file, *args):

        """ Write tandem repeat to file.

        Write tandem repeat to file using one of two formats.

        Args:
            format (str):  Either "stockholm" or "fasta"
            file (str): Path to output file

        .. todo:: Write checks for ``format`` and ``file``.

        """


        if format == 'stockholm':
            repeat_io.save_repeat_stockholm(self.msa, file)
        elif format == 'fasta':
            repeat_io.save_repeat_fasta(self.msa, file)
        else:
            raise Exception('format is unknown.')

    def score(self, score_type):
        if not hasattr(self,'dScore') or not score_type in self.dScore:
            self.calculate_scores(scoreslist=[score_type])
        return self.dScore[score_type]

    def pValue(self, score_type):
        if not hasattr(self,'dPValue') or not score_type in self.dPValue:
            self.calculate_pValues(scoreslist=[score_type])
        return self.dPValue[score_type]

    def divergence(self, score_type):
        if not hasattr(self,'dDivergence') or not score_type in self.dDivergence:
            self.calculate_scores(scoreslist=[score_type])
        return self.dDivergence[score_type]

    def __init__(self, msa, begin = None,
                 sequence_type = config["sequence_type"],
                 scoreslist=config['scoreslist'],
                 calc_score = config.as_bool('calc_score'), calc_pValue = config.as_bool('calc_pValue')):

        """
        .. todo:: if calc_score == False and calc_pValue == True, there is an error.
            However, calc_pValue == True should automatically require the score to be
            calculated.

        """

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
                self.calculate_scores(scoreslist=scoreslist)

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

    def calculate_scores(self, scoreslist=config['scoreslist']):

        """ Calculate scores on a Repeat instances for all scores in `scoreslist`.

        Args:
            self (Repeat): An instance of the repeat class.

        Kwargs:
            scoreslist(list of str): A list of the names of the scores (e.g. ["phylo_gap01"])
                that are calculated on the Repeat instance `self`
        """

        if not hasattr(self,'dScore'):
            self.dScore = defaultdict(float)

        if not hasattr(self,'dDivergence'):
            self.dDivergence = defaultdict(float)

        worst_score = {
            'entropy' : 1,
            'parsimony' : 1,
            'pSim' : 0,
            'phylo' : 0,
            'phylo_gap01': 0,
            'phylo_gap001': 0}

        if  self.lD == 0: # Enter worst score
            self.dScore = {iScore: worst_score[iScore] for iScore in scoreslist}
        elif self.n < 2: # Enter None, as a score cannot be calculated if there is just one repeat unit.
            self.dScore = {iScore: None for iScore in scoreslist}
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.dScore['entropy'] = repeat_score.meanSimilarity(self, repeat_score.entropy)
            if 'parsimony' in scoreslist:
                self.dScore['parsimony'] = np.round_(repeat_score.meanSimilarity(self, repeat_score.parsimony), decimals = config.as_int('precision'))
            if 'pSim' in scoreslist:
                self.dScore['pSim'] = np.round_(repeat_score.meanSimilarity(self, repeat_score.pSim), decimals = config.as_int('precision'))
            if 'phylo' in scoreslist:
                self.dDivergence['phylo'], self.dScore['phylo'] = repeat_score.phyloStarTopology_local(self)
            if 'phylo_gap01' in scoreslist:
                self.dDivergence['phylo_gap01'], self.dScore['phylo_gap01'] = repeat_score.phyloStarTopology_local(self, gaps = 'row_wise', indelRatePerSite = 0.01)
                self.dDivergence['phylo_gap01_ignore_trailing_gaps'], self.dScore['phylo_gap01_ignore_trailing_gaps'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps', indelRatePerSite = 0.01)
                self.dDivergence['phylo_gap01_ignore_coherent_deletions'], self.dScore['phylo_gap01_ignore_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_coherent_deletions', indelRatePerSite = 0.01)
                self.dDivergence['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'], self.dScore['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps_and_coherent_deletions', indelRatePerSite = 0.01)
            if 'phylo_gap001' in scoreslist:
#                self.dDivergence['phylo_gap001'], self.dScore['phylo_gap001'] = repeat_score.phyloStarTopology_local(self, gaps = 'row_wise', indelRatePerSite = 0.001)
#                self.dDivergence['phylo_gap001_ignore_trailing_gaps'], self.dScore['phylo_gap001_ignore_trailing_gaps'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps', indelRatePerSite = 0.001)
#                self.dDivergence['phylo_gap001_ignore_coherent_deletions'], self.dScore['phylo_gap001_ignore_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_coherent_deletions', indelRatePerSite = 0.001)
                self.dDivergence['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'], self.dScore['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = repeat_score.phyloStarTopology_local(self, gaps = 'ignore_trailing_gaps_and_coherent_deletions', indelRatePerSite = 0.001)

    def calculate_pValues(self, scoreslist=config['scoreslist']):

        """ Calculate pValues on a Repeat instances for all scores in `scoreslist`.

        Args:
            self (Repeat): An instance of the repeat class.

        Kwargs:
            scoreslist(list of Str): A list of the names of the scores (e.g. ["phylo_gap01"])
                for which p-Values are calculated on the Repeat instance `self`
        """

        if not hasattr(self,'dPValue'):
            self.dPValue = defaultdict(int)

        if  self.lD == 0: # Enter worst p-value
            self.dPValue = {iScore: 1.0 for iScore in scoreslist}
        elif self.n < 2: # Enter None, as a pValue cannot be calculated if there is just one repeat unit.
            self.dPValue = {iScore: None for iScore in scoreslist}
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.dPValue['entropy'] = repeat_pvalue.pValueFromEmpiricialList(self, 'entropy')
            if 'parsimony' in scoreslist:
                self.dPValue['parsimony'] = repeat_pvalue.pValuePars(self)
            if 'pSim' in scoreslist:
                self.dPValue['pSim'] = repeat_pvalue.pValuePSim(self)
            if 'phylo' in scoreslist:
                self.dPValue['phylo'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo')
            if 'phylo_gap' in scoreslist:
                self.dPValue['phylo_gap'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap')
            if 'phylo_gap01' in scoreslist:
                self.dPValue['phylo_gap01'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01')
                self.dPValue['phylo_gap01_ignore_trailing_gaps'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score('phylo_gap01_ignore_trailing_gaps'))
                self.dPValue['phylo_gap01_ignore_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score('phylo_gap01_ignore_coherent_deletions'))
                self.dPValue['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap01', self.score('phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'))
            if 'phylo_gap001' in scoreslist:
#                self.dPValue['phylo_gap001'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001')
#                self.dPValue['phylo_gap001_ignore_trailing_gaps'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score('phylo_gap001_ignore_trailing_gaps'))
#                self.dPValue['phylo_gap001_ignore_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score('phylo_gap001_ignore_coherent_deletions'))
                self.dPValue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = repeat_pvalue.pValueFromEmpiricialList(self, 'phylo_gap001', self.score('phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'))


    def deleteInsertionColumns(self):

        """ Create the tandem repeat attributes `msa`, `masT`, `l` and `n` without
            insertion columns.

        We define insertion columns as columns where there are more or equal gaps
        compared to chars.
        These columns are removed from both `msa` to create `msaD` and `msaT` to create
        `msaTD`.

        todo:: Check: is `totD` needed anywhere else? Is `textD` needed?
        """

        ## We define insertion columns as columns where there are more or equal gaps
        ## compared to chars.
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

        """ Calculate for each column in the transposed multiple repeat unit alignment
            `msaTD`

            *   how many deletions of what length begin in this column?:
                 `deletion_lengths`[`iColumn`] = list of deletion lengths
            *   how many insertions of what length begin in potential insertion columns
                 before this column?: `insertion_lengths`[`iColumn`] = list of insertion lengths

        Args:
            self (Repeat): An instance of the repeat class.
        """

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

    def repeat_in_sequence(self,sequence):

        """ Sanity check whether the `repeat` is part of the `sequence` (in which it
            was detected). In case, calculate the position of the `repeat` within the `sequence`.

            If yes: Return True, set self.begin to corrected value if necessary.
            If no: Return False.

        Args:
            sequence (sequence): A sequence instance.

        todo:: None-straight forward replaces such as "U" -> "C" are implemented. These
            need generalisation.
        todo:: save_original_msa is needed here?
        """

        repeat_sequence = "".join(self.msa).upper().replace("_", "").replace("-", "")
        starts = [m.start()+1 for m in re.finditer(repeat_sequence,sequence.replace("U", "C").replace("O", "K"))] # The first letter in the sequence is counted as 1 (not 0, as in python).

        if len(starts) != 0: # Is the tandem repeat predicted correctly?
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

########################################### MAIN #########################################

if __name__=="__main__":

    myTR = Repeat(begin = 10, msa = ['DFG','DFG','DFG'], sequence_type = 'AA', calc_score = True, calc_pValue = True)


