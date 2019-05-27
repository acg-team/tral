# (C) 2011, Alexander Korsunsky
# (C) 2011-2015 Elke Schaper

"""
    :synopsis: The Repeat Class.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

from collections import defaultdict
from copy import deepcopy
import logging
import numpy as np
import pickle
import re

from tral.repeat import repeat_score, repeat_pvalue, repeat_io, repeat_align
from tral import configuration

LOG = logging.getLogger(__name__)

CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat"]

# ################ Repeat class ###############################################


class Repeat:

    """ Tandem repeats are consecutive copies of genomic sequence.

    A tandem repeat is stored in terms of a multiple sequence alignment of its
    tandem repeat units. For example, the amino acid tandem repeat
    KLMKLMKLKLM
    can be displayed as a multiple sequence alignment (here: "multiple repeat
    unit alignment")
    KLM
    KLM
    KL-
    KLM
    and can be stored as
    ["KLM", "KLM", "KL-", "KLM"].
    The - char stands for a gap: It indicates that either, there used to be a
    char at the same position, but it got deleted. Or, all other chars in the
    same column got inserted at some point.

    For research on tandem repeats, many different statistics need to be
    calculated on the multiple sequence alignment of the tandem repeat.

    [Description of these statistics]

    Attributes:
        msa_original (list of str): The original alignment of repeat units.
        msa (list of str):  The alignment of repeat units without rare chars.
        msaT (list of str): The transposed alignment of repeat units. E.g.
            ["KKKK", "LLLL", "MM-M"].
        sequence_type (str): "AA" or "DNA"
        n (int): The number of repeat units in the original msa. E.g. 4.
        l_msa (int): The length of the repeat units in the original msa.
                     E.g. 3.
        l_effective (int): The number of columns in the original msa not
                           counting insertion columns (i.e. columns with more
                           gaps than chars).
        n_effective (float): The number of chars in the MSA without insertion
                             columns divided by n_effective.
        text (str): The repeat region of the tandem repeats (e.g. all chars in
                    the MSA without gaps concatenated.)
        n_gaps (int): The number of gaps in ``msa``
        repeat_region_length (int): The lengths of *text*
    """

    def __str__(self):
        """
            Create string for Repeat instance.
        """
        first_line = ">"
        if hasattr(self, 'begin'):
            first_line += " begin:{0}".format(self.begin)
        if hasattr(self, 'l_effective'):
            first_line += " l_effective:{0}".format(self.l_effective)
        if hasattr(self, 'n'):
            first_line += " n:{0}".format(self.n)
        if hasattr(self, 'd_pvalue') and len(self.d_pvalue) >= 1:
            if CONFIG['scoreslist'][-1] in self.d_pvalue:
                pvalue = self.d_pvalue[CONFIG['scoreslist'][-1]]
                model = CONFIG['scoreslist'][-1]
            else:
                model, pvalue = next(iter(self.d_pvalue.items()))
            try:
                divergence = self.d_divergence[model]
            except:
                divergence = None
            first_line += " pvalue:{0} divergence:{1}" \
                          " type:{2}".format(pvalue, divergence, model)

        return first_line + "\n" + "\n".join(self.msa)

    def create(file, input_format, **kwargs):
        """ Read tandem repeat from file.

        Read tandem repeat from file (currently, only pickle is supported)

        Args:
            format (str):  Currently only "pickle"
            file (str): Path to output file

        .. todo:: Write checks for ``format`` and ``file``.
        .. todo:: Write tests.
        .. todo:: Can repeat_io.read_fasta be replaced by e.g. BioPython seqIO?
        """

        if input_format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        elif input_format == 'fasta':
            l_repeat = repeat_io.read_fasta(file, **kwargs)
            return [Repeat(msa=i_msa, sequence_type=i_sequence_type)
                    for i_msa, i_sequence_type in l_repeat]
        elif input_format == 'simulate':
            l_repeat = repeat_io.random_sequence(return_type='repeat', **kwargs)
            return [Repeat(msa=i_msa, sequence_type=i_sequence_type)
                    for i_msa, i_sequence_type in l_repeat]
        else:
            raise Exception('input_format is unknown: {}.'.format(input_format))

    def write(self, file, file_format):
        """ Write tandem repeat to file.

        Write tandem repeat to file using one of three formats.

        Args:
            file (str): Path to input file
            file_format (str):  Either "fasta", "pickle" or "stockholm"

        .. todo:: Write checks for ``format`` and ``file``.

        """

        if file_format == 'fasta':
            repeat_io.save_repeat_fasta(self.msa, file)
        elif file_format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        elif file_format == 'stockholm':
            repeat_io.save_repeat_stockholm(self.msa, file)
        else:
            raise Exception('file_format is unknown: {}.'.format(file_format))

    def score(self, score_type):
        if not hasattr(self, 'd_score') or score_type not in self.d_score:
            self.calculate_scores(scoreslist=[score_type])
        if score_type in self.d_score:
            return self.d_score[score_type]

    def pvalue(self, score_type):
        if not hasattr(self, 'd_pvalue') or score_type not in self.d_pvalue:
            self.calculate_pvalues(scoreslist=[score_type])
        if score_type in self.d_pvalue:
            return self.d_pvalue[score_type]

    def divergence(self, score_type):
        if not hasattr(
                self,
                'd_divergence') or score_type not in self.d_divergence:
            self.calculate_scores(scoreslist=[score_type])
        if score_type in self.d_divergence:
            return self.d_divergence[score_type]

    def __init__(self, msa, begin=None,
                 sequence_type=CONFIG_GENERAL["sequence_type"],
                 scoreslist=CONFIG['scoreslist'],
                 calc_score=CONFIG.as_bool('calc_score'),
                 calc_pvalue=CONFIG.as_bool('calc_pvalue')):
        """
        .. todo:: if calc_score == False and calc_pvalue == True, there is an
            error. However, calc_pvalue == True should automatically require
            the score to be calculated.

        """

        LOG.debug(msa)

        # The index of begin start out on Zero.
        if begin is not None:
            self.begin = begin

        # Replace underscores and asterisks with scores, transform letters to
        # uppercase, replace Selenocysteine (U) with Cysteine (C),
        # replace Pyrrolysine (O) with Lysine (K).
        self.msa = [repeat_unit.upper().replace("_", "-").replace("*", "-")
                    for repeat_unit in msa]

        self.sequence_type = sequence_type

        self.msa_standard_aa = [standardize(i, self.sequence_type) for i in self.msa]

        self.n = len(self.msa)
        self.l_msa = max(map(len, self.msa))

        # Assure every line is equally long
        for i, iMSA in enumerate(self.msa):
            self.msa[i] += "-" * (self.l_msa - len(iMSA))
        LOG.debug(
            "Constructing repeat from MSA:\n%s",
            "\n".join(self.msa))

        # Transpose MSA
        self.msaT = ["".join(c) for c in zip(*self.msa)]

        # Get text for MSA.
        self.text = " ".join(self.msa)

        self.n_gaps = self.text.count("-")
        self.repeat_region_length = len(
            self.text) - self.n_gaps - self.text.count(" ")

        self.delete_insertion_columns()
        if self.l_effective != 0:
            self.gap_structure()

            if calc_score:
                self.calculate_scores(scoreslist=scoreslist)

            if calc_pvalue:
                self.calculate_pvalues(scoreslist=scoreslist)

    def calc_n_effective(self):
        """ Calculate the effective number of repeat units n_effective

        Calculate the number of letters in in all non-insertion columns,
        divided by <l_effective>

        Args:
            self (Repeat instance)

        """

        if self.l_effective != 0:
            self.n_effective = (
                len("".join(self.msaD)) - self.textD_standard_aa.count('-')) / self.l_effective
        else:
            self.n_effective = 0

    def calc_index_msa(self):
        ''' formerly named: set_begin'''

        # Calculate the transposed MSA with the index in protein sequence
        # instead of letter to be able to compare two overlapping Repeats
        # within the same sequence region.
        # python 2: range returns list by default
        index = list(
            range(
                self.begin + self.repeat_region_length - 1, self.begin - 1, -1))
        # Replace every letter with the index in the protein sequence
        msaI = [[j if j == '-' else index.pop() for j in i]
                for i in self.msa]  # the 'I' stands for index
        # Transpose msaI using list-comprehensions/zip/join magic
        self.msaIT = [[i for i in c if i != '-'] for c in zip(*msaI)]
        self.msaIT = [i for i in self.msaIT if len(i) > 0]

    def calculate_scores(self, scoreslist=CONFIG['scoreslist']):
        """ Calculate scores on a Repeat instance.

        Calculate scores on a Repeat instance for all scores in `scoreslist`.

        Args:
            self (Repeat): An instance of the repeat class.

        Kwargs:
            scoreslist(list of str): A list of the names of the scores (e.g.
             ["phylo_gap01"]) that are calculated on the Repeat instance `self`
        """

        if not hasattr(self, 'd_score'):
            self.d_score = defaultdict(float)

        if not hasattr(self, 'd_divergence'):
            self.d_divergence = defaultdict(float)

        worst_score = {
            'entropy': 1,
            'parsimony': 1,
            'pSim': 0,
            'phylo': 0,
            'phylo_gap01': 0,
            'phylo_gap001': 0}

        if self.l_effective == 0:  # Enter worst score
            self.d_score = {
                iScore: worst_score[iScore] for iScore in scoreslist}
        # Enter None, as a score cannot be calculated if there is just one
        # repeat unit.
        elif self.n < 2:
            self.d_score = {iScore: None for iScore in scoreslist}
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.d_score['entropy'] = \
                    repeat_score.mean_similarity(self, repeat_score.entropy)
            if 'parsimony' in scoreslist:
                self.d_score['parsimony'] = np.round(
                    repeat_score.mean_similarity(
                        self,
                        repeat_score.parsimony),
                    decimals=CONFIG.as_int('precision'))
            if 'pSim' in scoreslist:
                self.d_score['pSim'] = np.round(
                    repeat_score.mean_similarity(
                        self,
                        repeat_score.pSim),
                    decimals=CONFIG.as_int('precision'))
            if 'phylo' in scoreslist:
                self.d_divergence['phylo'], self.d_score['phylo'] = \
                    repeat_score.phylo_star_topology_local(self)
            if 'phylo_gap01' in scoreslist:
                self.d_divergence['phylo_gap01'], self.d_score['phylo_gap01'] = \
                    repeat_score.phylo_star_topology_local(self,
                                                           gaps='row_wise',
                                                           indel_rate_per_site=0.01)
            if 'phylo_gap01_ignore_trailing_gaps' in scoreslist:
                self.d_divergence['phylo_gap01_ignore_trailing_gaps'], \
                    self.d_score['phylo_gap01_ignore_trailing_gaps'] = \
                    repeat_score.phylo_star_topology_local(self, gaps='ignore_trailing_gaps', indel_rate_per_site=0.01)
            if 'phylo_gap01_ignore_coherent_deletions' in scoreslist:
                self.d_divergence['phylo_gap01_ignore_coherent_deletions'], self.d_score[
                    'phylo_gap01_ignore_coherent_deletions'] = repeat_score.phylo_star_topology_local(self, gaps='ignore_coherent_deletions', indel_rate_per_site=0.01)
            if 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions' in scoreslist:
                self.d_divergence['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'], \
                    self.d_score['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = \
                    repeat_score.phylo_star_topology_local(
                        self, gaps='ignore_trailing_gaps_and_coherent_deletions', indel_rate_per_site=0.01)
            if 'phylo_gap001' in scoreslist:
                self.d_divergence['phylo_gap001'], self.d_score['phylo_gap001'] = repeat_score.phylo_star_topology_local(
                    self, gaps='row_wise', indel_rate_per_site=0.001)
            if 'phylo_gap001_ignore_trailing_gaps' in scoreslist:
                self.d_divergence['phylo_gap001_ignore_trailing_gaps'], self.d_score[
                    'phylo_gap001_ignore_trailing_gaps'] = repeat_score.phylo_star_topology_local(self, gaps='ignore_trailing_gaps', indel_rate_per_site=0.001)
            if 'phylo_gap001_ignore_coherent_deletions' in scoreslist:
                self.d_divergence['phylo_gap001_ignore_coherent_deletions'], self.d_score[
                    'phylo_gap001_ignore_coherent_deletions'] = repeat_score.phylo_star_topology_local(self, gaps='ignore_coherent_deletions', indel_rate_per_site=0.001)
            if 'phylo_gap001_ignore_trailing_gaps_and_coherent_deletions' in scoreslist:
                self.d_divergence['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'], self.d_score[
                    'phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = repeat_score.phylo_star_topology_local(self, gaps='ignore_trailing_gaps_and_coherent_deletions', indel_rate_per_site=0.001)

    def calculate_pvalues(self, scoreslist=CONFIG['scoreslist']):
        """ Calculate pvalues on a Repeat instance.

        Calculate pvalues on a Repeat instance for all scores in `scoreslist`.

        Args:
            self (Repeat): An instance of the repeat class.

        Kwargs:
            scoreslist(list of Str): A list of the names of the scores (e.g.
            ["phylo_gap01"]) for which p-Values are calculated on the
            Repeat instance `self`
        """
        if type(scoreslist) != list:
            raise ValueError("Please make sure that scoreslist is a list. " \
                            "You can make this sure by a comma at the end of your Config-File.")
        
        if not hasattr(self, 'd_pvalue'):
            self.d_pvalue = defaultdict(int)

        if self.l_effective == 0:  # Enter worst p-value
            self.d_pvalue = {iScore: 1.0 for iScore in scoreslist}
        # Enter None, as a pvalue cannot be calculated if there is just one
        # repeat unit.
        elif self.n < 2:
            self.d_pvalue = {iScore: None for iScore in scoreslist}
        else:   # Calculate score
            if 'entropy' in scoreslist:
                self.d_pvalue['entropy'] = \
                    repeat_pvalue.pvalue_from_empirical_list(self, 'entropy')
            if 'parsimony' in scoreslist:
                self.d_pvalue['parsimony'] = repeat_pvalue.pvalue_pars(self)
            if 'pSim' in scoreslist:
                self.d_pvalue['pSim'] = repeat_pvalue.pvalue_psim(self)
            if 'phylo' in scoreslist:
                self.d_pvalue['phylo'] = \
                    repeat_pvalue.pvalue_from_empirical_list(self, 'phylo')
            if 'phylo_gap' in scoreslist:
                self.d_pvalue['phylo_gap'] = \
                    repeat_pvalue.pvalue_from_empirical_list(self, 'phylo_gap')
            if 'phylo_gap01' in scoreslist:
                self.d_pvalue['phylo_gap01'] = \
                    repeat_pvalue.pvalue_from_empirical_list(self, 'phylo_gap01')
            if 'phylo_gap01_ignore_trailing_gaps' in scoreslist:
                self.d_pvalue['phylo_gap01_ignore_trailing_gaps'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                        self,
                        'phylo_gap01',
                        self.score('phylo_gap01_ignore_trailing_gaps'))
            if 'phylo_gap01_ignore_coherent_deletions' in scoreslist:
                self.d_pvalue['phylo_gap01_ignore_coherent_deletions'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                        self,
                        'phylo_gap01',
                        self.score('phylo_gap01_ignore_coherent_deletions'))
            if 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions' in scoreslist:
                self.d_pvalue['phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                        self,
                        'phylo_gap01',
                        self.score('phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'))
            if 'phylo_gap001' in scoreslist:
                self.d_pvalue['phylo_gap001'] = \
                    repeat_pvalue.pvalue_from_empirical_list(self,
                                                             'phylo_gap001')
            if 'phylo_gap001_ignore_trailing_gaps' in scoreslist:
                self.d_pvalue['phylo_gap001_ignore_trailing_gaps'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                    self,
                    'phylo_gap001',
                    self.score('phylo_gap001_ignore_trailing_gaps'))
            if 'phylo_gap001_ignore_coherent_deletions' in scoreslist:
                self.d_pvalue['phylo_gap001_ignore_coherent_deletions'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                    self,
                    'phylo_gap001',
                    self.score('phylo_gap001_ignore_coherent_deletions'))
            if 'phylo_gap001_ignore_trailing_gaps_and_coherent_deletions' in scoreslist:
                self.d_pvalue['phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'] = \
                    repeat_pvalue.pvalue_from_empirical_list(
                    self,
                    'phylo_gap001',
                    self.score('phylo_gap001_ignore_trailing_gaps_and_coherent_deletions'))

    def delete_insertion_columns(self):
        """ Create the tandem repeat attributes `msa`, `masT`, `l` and `n`
            without insertion columns.

        We define insertion columns as columns where there are more or equal
        gaps compared to chars.
        These columns are removed from both `msa` to create `msaD` and `msaT`
        to create `msaTD`.

        """

        # We define insertion columns as columns where there are more or equal
        # gaps compared to chars.
        self.msaTD = [column for column in self.msaT
                      if (2 * column.count('-') < len(column))
                      ]
        self.insertions = re.compile(r"-+").findall(''.join([
            'p' if (2 * column.count('-') < len(column)) else '-'
            for column in self.msaT
        ])
        )

        self.msaTD_standard_aa = [standardize(i, self.sequence_type) for i in self.msaTD]

        self.msaD = ["".join(c) for c in zip(*self.msaTD)]

        self.l_effective = len(self.msaTD)
        textD = " ".join(self.msaD)
        self.textD_standard_aa = standardize(textD, self.sequence_type)
        # totD is used in repeat_score
        self.totD = len(textD) - textD.count('-') - self.n + 1
        if self.l_effective != 0:
            self.n_effective = (
                len("".join(self.msaD)) - textD.count('-')) / self.l_effective
        else:
            self.n_effective = 0

    def gap_structure(self):
        ''' Calculate the number and length of insertions and deletions for
            this tandem repeat.'''

        # Needed: Number and size of insertions, number and size of deletions
        # row_wise: (standard case) column wise coherence is ignored, the
        # indels are summarised over the columns.
        # ignore_trailing_gaps: In addition to the standard case, trailing gaps
        # on either side of the tandem repeat are replaced with letters, and
        # hence ignored.
        # The separation into insertions and deletions comes before ignoring
        # trailing gaps.
        # ignore_coherent_deletions : In addition to the standard case, same
        # size deletions spanning the same columns are only counted once.
        self.insertions = []
        self.deletions = {}
        self.gaps = {}

        # row_wise
        # 1. Detect insertions
        # First, name columns either '-' if they are insertion columns, or 'p',
        #  if they are not.
        # Concetanate the string of 'p's and '-', and double it, to consider
        # insertions in the last column
        # and the first column as associated (= of the same origin)
        # Lastly, detect the length of all insertions.
        insertions = re.compile(r"-+").findall(2 * (''.join([
            'p' if (2 * column.count('-') < len(column)) else '-'
            for column in self.msaT
        ])))
        insertions = insertions[int(len(insertions) / 2):-1] if len(
            insertions) % 2 == 1 else insertions[int(len(insertions) / 2):]
        self.insertions = [len(i) for i in insertions]
        if self.l_effective == 0:
            self.insertions = [self.l_msa]

        # 2. Detect deletions
        # CHECK this lines. you used msaTD before, but swapped to msaD
        deletions = [(m.start() % self.l_msa, len(m.group()))
                     for m in re.finditer(re.compile(r"-+"), "".join(self.msaD))]
        self.deletions['row_wise'] = [i[1] for i in deletions]
        self.gaps['row_wise'] = self.insertions + self.deletions['row_wise']

        # ignore_coherent_deletions
        # Remove duplicates from the list of tuples(index of deletion, deletion
        # length)
        deletions_ignore_coherent = {iD: None for iD in deletions}
        self.deletions['ignore_coherent_deletions'] = [i[1]
                                                       for i in deletions_ignore_coherent.keys()]
        self.gaps['ignore_coherent_deletions'] = self.insertions + \
            self.deletions['ignore_coherent_deletions']

        # ignore_trailing_gaps
        # To be precise, we should ignore_trailing_gaps only, if we have
        # ascertained that there were no letters in possible insertion columns
        # before the first deletion column. For now, this case is ignored.
        msaD = deepcopy(self.msaD)

        def sub_trailing_gap(match_object):
            return 'x' * len(match_object.group())
        msaD[0] = re.sub('^-+', sub_trailing_gap, msaD[0])
        msaD[-1] = re.sub('-+$', sub_trailing_gap, msaD[-1])

        # CHECK this lines. you used msaTD before, but swapped to msaD
        # msaTD = ["".join(c) for c in zip(*msaD)]
        deletions = [(m.start() % self.l_msa, len(m.group()))
                     for m in re.finditer(re.compile(r"-+"), "".join(msaD))]
        self.deletions['ignore_trailing_gaps'] = [i[1] for i in deletions]
        self.gaps['ignore_trailing_gaps'] = self.insertions + \
            self.deletions['ignore_trailing_gaps']

        ## ignore_trailing_gaps & ignore_coherent_deletions
        deletions_ignore_coherent = {iD: None for iD in deletions}
        self.deletions['ignore_trailing_gaps_and_coherent_deletions'] = [
            i[1] for i in deletions_ignore_coherent.keys()]
        self.gaps['ignore_trailing_gaps_and_coherent_deletions'] = self.insertions + \
            self.deletions['ignore_trailing_gaps_and_coherent_deletions']

    def gap_structure_HMM(self):
        """ Calculate for each column in the transposed multiple repeat unit
            alignment `msaTD`

            *   how many deletions of what length begin in this column?:
                 `deletion_lengths`[`iColumn`] = list of deletion lengths
            *   how many insertions of what length begin in potential insertion
                columns before this column?:
                `insertion_lengths`[`iColumn`] = list of insertion lengths

        Args:
            self (Repeat): An instance of the repeat class.
        """

        # Insertions: In which column (with reference to self.msaTD) do how
        # many insertions of what length start?
        # In which columns are insertions?
        insertions = [
            i for i,
            column in enumerate(
                self.msaT) if (
                2 *
                column.count('-') >= len(column))]

        # Which insertions are neighboring == form a block, and on which index
        # with reference to self.msaTD do they start?
        insertion_blocks = {}
        n_insertion_columns = 0
        while insertions:
            current = insertions.pop(0)
            insertion_blocks[current] = {'indices_self.msaT': [current]}
            while insertions and insertions[0] == insertion_blocks[
                    current]['indices_self.msaT'][-1] + 1:
                insertion_blocks[current][
                    'indices_self.msaT'].append(insertions.pop(0))
            insertion_blocks[current]['index_self.msaTD'] = insertion_blocks[
                current]['indices_self.msaT'][0] - n_insertion_columns
            n_insertion_columns += len(
                insertion_blocks[current]['indices_self.msaT'])

        # If an insertions ranges over the repeat unit border, take the
        # insertion on both parts together:
        if 0 in insertion_blocks.keys() and len(
                self.msaT) - 1 in insertion_blocks.keys():
            insertion_blocks[current]['indices_self.msaT'].extend(
                insertion_blocks[0]['indices_self.msaT'])
            insertion_blocks[current][
                'index_self.msaTD'] = insertion_blocks[0]['index_self.msaTD']
            del insertion_blocks[0]

        # Restructure insertion_blocks: We are only interested in the index
        # with reference to self.msaTD, and the start and end index with
        # reference to self.msaT:
        insertion_blocks = {
            i['index_self.msaTD']: (
                i['indices_self.msaT'][0],
                i['indices_self.msaT'][
                    -1]) for i in insertion_blocks.values()}

        # Go to each insertion block, and see how many insertion it includes,
        # and of what length.
        l = len(self.msaT)
        insertion_lengths = {}
        for iL in range(self.l_effective):
            insertion_lengths[iL] = []
            if iL in insertion_blocks:
                begin = insertion_blocks[iL][0]
                end = insertion_blocks[iL][1]
                length = end - begin + 1
                # In case the insertion blocks bridges over then end of the
                # repeat unit:
                if begin > end:
                    length += l
                print(length)

                repeat_sequence = get_repeat_sequence(self.msa)
                for i in range(self.n):
                    insertion = repeat_sequence[
                        i *
                        l +
                        begin:min(
                            len(repeat_sequence),
                            i *
                            l +
                            begin +
                            length)]
                    insertion_lengths[iL].extend(
                        [len(iI) for iI in re.findall(r'\w+', insertion)])

                # In case the insertion blocks bridges over then end of the
                # repeat unit:
                if begin > end:
                    insertion = repeat_sequence[:end + 1]
                    insertion_lengths[iL].extend(
                        [len(iI) for iI in re.findall(r'\w+', insertion)])

        LOG.debug(
            "The insertion_lengths of the current TR are: {0}".format(insertion_lengths))

        # Deletions are much easier to handle than insertions, as they are
        # already in the right coordinate system (self.msaTD)
        # Find all deletions, and note their start index, and their lengths.
        deletions_all = [(m.start(), len(m.group()))
                         for m in re.finditer('-+', ''.join(self.msaD))]
        # Group the deletions with respect to the column in which they start:
        deletion_lengths = {iL:
                            [iD[1] for iD in deletions_all if iD[0] % self.l_effective == iL]
                            for iL in range(self.l_effective)}
        LOG.debug(
            "The deletion_lengths of the current TR are: {0}".format(deletion_lengths))

        return insertion_lengths, deletion_lengths

    def save_original_msa(self, sequence):
        ''' recalculate the original msa, assuming that begin is defined
        correctly (index starting counting on 1.) This function might be
        obsolete when the original MSA is saved from the beginning.'''
        repeat_sequence = sequence[
            self.begin -
            1:self.begin -
            1 +
            self.repeat_region_length]

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
    
    def realign_TR(self, realignment='proPIP', rate_distribution=CONFIG['castor_parameter']['rate_distribution']):
        """ Realign multiple sequence alignment of tandem repeat

        Args:
            self (Repeat instance)
            realignment (str): either "proPIP" or "mafft"
            rate_distribution (str): either "constant" or "gamma" (per default value from REPEAT_CONFIG['castor_parameter']['rate_distribution'])
        """
        self.msa = repeat_align.realign_repeat(self.msa, 
                                            realignment, 
                                            sequence_type = self.sequence_type, 
                                            rate_distribution = rate_distribution)

# Standardize MSA #############################################################


def standardize(blob, sequence_type):

    for original, replacement in CONFIG_GENERAL[sequence_type][
            'ambiguous_chars'].items():
        blob = blob.replace(original, replacement[0])
    return blob


def get_repeat_sequence(msa):

    return "".join([i.replace("_", "").replace("-", "") for i in msa])


# Localize repeat #############################################################


def calculate_position_in_alignment(begin, length, alignment):
    """ Calculate the index of the begin and the end of a TR within an alignment
        returns

        Calculate the index of the begin and the end of a TR within an alignment
        returns

    Args:
        begin (int): (needs explaining)
        length (int): (needs explaining)
        alignment (?): (needs explaining)

    Returns:
        Dict
    """

    # a alignment.seq (hopefully!) is a string of letters, with many gaps "-",
    # as it this particular case is an alignment sequence.
    seq_indexed = [i for i, j in enumerate(str(alignment)) if j != '-']
    LOG.debug("begin: {0}; length: {1}".format(str(begin), str(length)))
    LOG.debug("alignment: {0}".format(str(alignment)))
    return {
        'begin': seq_indexed[
            begin -
            1],
        'end': seq_indexed[
            begin +
            length -
            2]}
