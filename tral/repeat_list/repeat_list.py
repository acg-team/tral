# (C) 2015 Elke Schaper

from tral.repeat_list import repeat_list_io
"""
    :synopsis: The RepeatList Class.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import pickle
import sys

LOG = logging.getLogger(__name__)


class RepeatList:

    """ A `RepeatList` contains list of repeats that belong to the same
    sequence, or set of sequences.

    `RepeatList` contains methods that act on several tandem repeats.
    For example, methods to

    *   detect overlapping tandem repeats
    *   identify the tandem repeat with highest statistical significance
    *   filter a set of tandem repeats according to different filtering
        procedures.


    Attributes:
        repeats (list of Repeat): The list of tandem repats.
    """

    def __str__(self):
        """
            Create string for RepeatList instance.
        """
        return str(self.repeats)

    def __init__(self, repeats):
        self.repeats = repeats

    def __add__(self, rl):
        if rl:
            return RepeatList(self.repeats + rl.repeats)
        else:
            return self

    def intersection(self, rl):
        if rl:
            return RepeatList([i for i in self.repeats if i in rl.repeats])
        else:
            return self

    def create(file, input_format):
        """ Read ``RepeatList`` from file.

        Read ``RepeatList`` from file (currently, only pickle is supported)

        Args:
            input_format (str):  Currently only "pickle"
            file (str): Path to output file

        .. todo:: Write checks for ``input_format`` and ``file``.

        """

        if input_format == 'pickle':
            with open(file, 'rb') as fh:
                return pickle.load(fh)
        else:
            raise Exception('input_format is unknown: {}.'.format(input_format))

    def write(self, output_format, file=None, return_string=None, *args):
        """ Serialize and write ``RepeatList`` instances.

        Serialize ``RepeatList`` instance using the stated ``output_format``.
        If a ``file`` is specified, save the String. If ``str`` is specified, give back
        the String (not possible for pickles).

        Args:
            output_format (str):  The output format: Either "pickle" or "tsv"
            file (str): Path to output file

        .. todo:: Write checks for ``output_format`` and ``file``.
        """

        if output_format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
            file = False
        elif output_format == 'tsv':
            output = repeat_list_io.serialize_repeat_list_tsv(self)
        else:
            raise Exception('output_format is unknown: {}'.format(output_format))

        if file:
            with open(file, 'w') as fh:
                fh.write(output)
        if return_string:
            return output

    def filter(self, func_name, *args, **kwargs):
        """ Filter ``repeats`` according to ``func_name``.

        Filter ``repeats`` according to ``func_name``

        Args:
            func_name (str): The name of a local filtering method.

        """

        # Check: is func_name a str, or is it a method? if it is a method, it
        # must be in dir(): ``funcname in dir`` == TRUE?

        func = getattr(sys.modules[__name__], func_name)
        return RepeatList(func(self, *args, **kwargs))

    def cluster(self, overlap_type, *args):
        """ Cluster ``repeats`` according to ``overlap_type``.

        Cluster ``repeats`` according to ``overlap_type``. We assume that overlap of
        repeats is transitive: If A overlaps with B, and B overlaps with C, (A,B,C) form
        one cluster.
        The attribute ``cluster`` is initiated to a dict:

            self.cluster = {"overlap_type1": [(0,2),(1)], "overlap_type2": [(0),(1), (2)]}

        In this toy example, the first and the third repeat in ``repeats`` cluster
        according to "overlap_type1", whereas no repeats cluster according to
        "overlap_type2".

        Args:
            overlap_type (str): The name of a local pairwise repeat overlapping method.

        """

        if not hasattr(self, 'd_cluster'):
            self.d_cluster = {}

        l_cluster = []

        l_repeat_indices = list(range(len(self.repeats)))
        while l_repeat_indices:
            i_cluster = [l_repeat_indices.pop()]
            l_repeat_in_cluster_check = i_cluster[:]

            while l_repeat_in_cluster_check:
                check = l_repeat_in_cluster_check.pop()
                remaining_indices = l_repeat_indices[:]
                for i in remaining_indices:
                    if two_repeats_overlap(
                            overlap_type,
                            self.repeats[check],
                            self.repeats[i]):
                        i_cluster.append(i)
                        l_repeat_in_cluster_check.append(i)
                        l_repeat_indices.remove(i)

            l_cluster.append(i_cluster)

        l_cluster = [set(i) for i in l_cluster]
        self.d_cluster[overlap_type] = l_cluster


# FILTERING
# Filter methods (func) are all defined in the same way: They return a subset
# of self.repeats, also as a list.
#  Return is a (list of `Repeat`).

def pvalue(rl, score, threshold):
    """ Returns all repeats in ``rl`` with a p-Value below a certain threshold.

    Returns all repeats in ``rl`` with a p-Value below a certain threshold.

    Args:
        rl (RepeatList): An instance of the RepeatList class.
        score (str): The type of score defines the pvalue that is used for filtering
        threshold (float): All repeats with a pvalue of type `score` above this threshold
            are filtered out.

    """

    threshold = float(threshold)

    res = []
    for i_repeat in rl.repeats:
        if i_repeat.pvalue(score) <= threshold:
            res.append(i_repeat)
    return res


def divergence(rl, score, threshold):
    """ Returns all repeats in ``rl`` with a divergence below a certain
        threshold.

    Returns all repeats in ``rl`` with a divergence below a certain threshold.

    Args:
        rl (RepeatList): An instance of the RepeatList class.
        score (str): The type of score defines the divergence that is used for
                    filtering
        threshold (float): All repeats with a divergence of type `score` above
                            this threshold are filtered out.
    """

    res = []
    for i_repeat in rl.repeats:
        if i_repeat.divergence(score) <= threshold:
            res.append(i_repeat)
    return res


def attribute(rl, attribute, type, threshold):
    """ Returns all repeats in ``rl`` with a attribute below (above) a certain
        threshold.

    Returns all repeats in ``rl`` with a attribute below (above) a certain
    threshold.

    Args:
        rl (RepeatList): An instance of the RepeatList class.
        attribute (str): The attribute of the Repeat instance
        type (str): Either "min" or "max"
        threshold (float): All repeats with an attribute value below (above)
                          this threshold are filtered out.
    """

    threshold = float(threshold)

    res = []
    for i_repeat in rl.repeats:
        value = getattr(i_repeat, attribute)
        if type == "min":
            if value >= threshold:
                res.append(i_repeat)
        elif type == "max":
            if value <= threshold:
                res.append(i_repeat)
        else:
            raise Exception(
                "type must either be 'min' or 'max'. Instead, it is {}.".format(type))
    return res


def none_overlapping_fixed_repeats(rl, rl_fixed, overlap_type):
    """ Returns all repeats in ``rl`` none-overlapping with ``rl_fixed``.

    Returns all repeats in ``rl`` none-overlapping with ``rl_fixed`` according
    to ``overlap``.

    Args:
        rl (RepeatList): An instance of the RepeatList class.
        rl (rl_fixed): A second instance of the RepeatList class.
        overlap (list): First list element: Name (str) of an overlap method in
                        repeat_list. All remaining elements are additional
                        arguments for this class.
    """

    res = []
    for i_repeat in rl.repeats:

        for i_repeat_fixed in rl_fixed.repeats:
            if two_repeats_overlap(
                    overlap_type,
                    repeat1=i_repeat,
                    repeat2=i_repeat_fixed):
                break
        else:
            res.append(i_repeat)

    return res


def none_overlapping(rl, overlap, l_criterion):
    """ Returns all none-overlapping repeats in ``rl``.

    Returns all none-overlapping repeats in ``rl``. Repeats are clustered
    according to ``overlap``. Of each cluster, only the best repeat is returned
    according to ``dCriterion``.

    Args:
        rl (RepeatList): An instance of the RepeatList class.
        overlap (tuple): First element: Name (str) of an overlap method in
            `repeat_list`. Second element: **kwargs. All remaining elements are
            additional arguments for this class.
        l_criterion (list): list of (criterion (str), criterion arguments)
            tuples. Until only one repeat is remainining in a cluster, the
            criteria are applied in order.
    """

    overlap_type = overlap[0]
    if len(overlap) > 1:  # has to be adapted to a more convenient way
        overlap_args = overlap[1]

    if not (hasattr(rl, 'd_cluster') and overlap_type in rl.d_cluster):
        if len(overlap) > 1:
            rl.cluster(overlap_type, overlap_args)
        else:
            rl.cluster(overlap_type)

    res = []
    for i_cluster in rl.d_cluster[overlap_type]:

        i_repeat = [rl.repeats[i] for i in i_cluster]

        for ctype, cvalue in l_criterion:
            criterion_type, criterion_value = ctype, cvalue
            if len(i_repeat) == 1:
                res.append(i_repeat[0])
                break

            if criterion_type == 'pvalue':
                min_value = min(i.pvalue(criterion_value) for i in i_repeat)
                i_repeat = [
                    i for i in i_repeat if i.pvalue(criterion_value) == min_value]

            elif criterion_type == 'divergence':
                min_value = min(i.divergence(criterion_value) for i in i_repeat)
                i_repeat = [
                    i for i in i_repeat if i.divergence(criterion_value) == min_value]

        else:
            LOG.debug(
                "repeat_list.none_overlapping(): > 1 Repeats have the same values...")
            res.append(i_repeat[0])

    return res


# OVERLAP DETECTION METHODS
def two_repeats_overlap(overlap_type, repeat1, repeat2):
    """ Helper method to test the overlap of ``repeat1`` and ``repeat2``.

    Helper method to test the overlap of ``repeat1`` and ``repeat2``. The
    overlap is calculated by the local method ``overlap_type``.

    Args:
        overlap_type (str): The name of a local pairwise repeat overlapping
                            method.
        repeat1 (Repeat): An instance of the Repeat class
        repeat2 (Repeat): A second instance of the Repeat class

    Returns:
        Forwards method output
    """

    is_overlapping = getattr(sys.modules[__name__], overlap_type)
    # Alternative implementation:
    # is_overlapping = globals()[overlap_type]
    return is_overlapping(repeat1, repeat2)


def shared_char(repeat1, repeat2):
    """ Do two TRs share at least one char?

    Return 1 if the two TRs share at least one char (amino acids or nucleotides); else 0.

    Args:
        repeat1 (Repeat): An instance of the Repeat class
        repeat2 (Repeat): A second instance of the Repeat class

    Returns:
        Bool: 1 if the repeats share >= 1 char, else 0.
    """

    if (repeat1.begin + repeat1.repeat_region_length - 1 < repeat2.begin) or \
       (repeat2.begin + repeat2.repeat_region_length - 1 < repeat1.begin):
        return False
    else:
        return True


def common_ancestry(repeat1, repeat2):
    """ Do two TRs share at least one pair of chars with common ancestry?

    Return 1 if the two TRs share at least one pair of chars (amino acids or
    nucleotides) with common ancestry; else 0.

    Args:
        repeat1 (Repeat): An instance of the Repeat class
        repeat2 (Repeat): A second instance of the Repeat class

    Returns:
        Bool: 1 if the repeats share >= 1 pair of chars with common ancestry, else 0.

    """

    if not hasattr(repeat1, 'msaIT'):
        repeat1.calc_index_msa()
    if not hasattr(repeat2, 'msaIT'):
        repeat2.calc_index_msa()

    original = repeat2.msaIT
    potential = repeat1.msaIT
    try:
        for p in potential:
            iP = 0
            for iP in range(len(p) - 1):
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
                        if p[iPRest] in original[i][j + 1:]:
                            #coverage = repeat1.repeat_region_length/repeat2.repeat_region_length
                            #greediness = repeat1.l_effective/repeat2.l_effective
                            return True
    except:
        LOG.warning('error in shared_char with original %s and potential %s',
                    str(repeat1.msaIT), str(repeat2.msaIT))
        LOG.warning('original:', str(repeat2.msa))
        LOG.warning('original:', str(repeat2.begin))
        LOG.warning('potential:', str(repeat1.msa))
        LOG.warning('potential:', str(repeat1.begin))
        return False
    return False
