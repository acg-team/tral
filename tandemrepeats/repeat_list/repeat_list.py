class Repeat_list:

    """ A `repeat_list` is a list of repeats that belong to the same sequence, or set of
        sequences.

    Repeat_set contains methods that act on several tandem repeats.
    For example, methods to

    *   detect overlapping tandem repeats
    *   identify the tandem repeat with highest statistical significance
    *   filter a set of tandem repeats according to different filtering procedures.


    Attributes:
        repeats (list of Repeat): The list of tandem repats.
    """


    def __init__(self, repeats):
        self.repeats = repeats

    def filter(self, func, *args):
        return Repeat_set( func(self.repeats, *args) )


### FILTER METHODS
# Filter methods are all defined in the same way: Their first input argument is a list of
# repeats. They return a subset of these repeats, also as a list.
#   First input is a (list of `Repeat`).
#   Return is a (list of `Repeat`).


def pValue(lRepeat, score, threshold):

    """ Repeat_list filter methods: Repeats with a p-Value below a certain threshold are
        filtered out.

    __init__ takes HMM parameters (including the alphabet, emission probabilities
    and transition probabilities) as input, and assigns them to class attributes.

    Args:
        lRepeat (list of Repeat): The unfiltered list of repeats
        score (str): The type of score defines the pValue that is used for filtering
        threshold (float): All repeats with a pValue of type `score` above this threshold
            are filtered out.

    """

    res = []
    for iRepeat in lRepeat:
        if iRepeat.pValues[score] <= threshold:
            res.append(x)
    return(res)
