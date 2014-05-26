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
# Needs documentation
# Requirement:
#    First input is a (list of `Repeat`).
#   Return is a (list of `Repeat`).


def pValue(lRepeat, score, threshold):
    res = []
    for iRepeat in lRepeat:
        if iRepeat.pValues[score] < threshold:
            res.append(x)
    return(res)
