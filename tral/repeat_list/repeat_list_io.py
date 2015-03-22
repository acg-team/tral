# (C) 2015 Elke Schaper


"""

    :synopsis: Input/output for the repeat_list class

    .. moduleauthor:: Elke Schaper <elke@inf.ethz.ch>

"""
import logging
import os

from tral import configuration
from tral.paths import *

log = logging.getLogger(__name__)

c = configuration.Configuration.Instance()
config_general = c.config
config = config_general["repeat_list"]


def serialize_repeat_list_tsv(tandem_repeats, config=config, *args):
    ''' Serialize a ``repeat_list`` instance as tsv.

        config defines which tandem repeat characteristics are added to the .tsv.
        This could be for example:

        * begin: position of the tandem repeats within the sequence,
        * pvalue: statistical significance of the tandem repeats
        * divergence: divergence of the tandem repeat units
        * l_effective: length of the tandem repeat units
        * n_effective: number of tandem repeat units

    Attributes:
        tandem_repeats (repeat_list): A ``Repeat_list`` instance.
        config (dict): A dictionary. E.g.: {"output_characteristics": ["n_effective", "pvalue"], model = "phylo_gap01"}
        *args: Additional arguments

    Returns:
        str: The tsv as a string.

    '''

    output_characteristics = config["output_characteristics"]
    model = config["model"]

    data = ["\t".join(output_characteristics)]
    for iRepeat in tandem_repeats.repeats:
        d = []
        for iCharacteristic in output_characteristics:
            if "divergence" == iCharacteristic:
                if hasattr(iRepeat, "dScore") and model in iRepeat.dScore:
                    d.append(iRepeat.divergence(model))
                else:
                    d.append(None)
            elif "pvalue" == iCharacteristic:
                if hasattr(iRepeat, "dPValue") and model in iRepeat.dPValue:
                    d.append(iRepeat.pvalue(model))
                else:
                    d.append(None)
            elif "score" == iCharacteristic:
                if hasattr(iRepeat, "dScore") and model in iRepeat.dScore:
                    d.append(iRepeat.score(model))
                else:
                    d.append(None)
            elif "msa_original" == iCharacteristic:
                try:
                    d.append(",".join(iRepeat.msa_original))
                except:
                    raise Exception(
                        "The attribute msa_original is not available for iRepeat.")
            elif getattr(iRepeat, iCharacteristic):
                d.append(getattr(iRepeat, iCharacteristic))
            else:
                raise Exception(
                    "The attribute {} is not available for tandem_repeats".format(iCharacteristic))
        data.append("\t".join(str(i) for i in d))

    return "\n".join(data) + "\n"


def save_repeat_fasta(tandem_repeats, file):
    ''' save multiple <tandem_repeats> in Fasta format in specified <file>

        At current, only one TR per sequence can be defined, as the identifiers in
        the dict <tandem_repeats> must be unique.

        Attributes:
            tandem_repeats: Dict of tandem repeats and identifiers. e.g. {'ENSP00012': msa1, 'ENSP00013': msa2}

            >ID
            GHKI
            GHKI
            GH--
    '''

    with open(file, 'w', newline='\n') as f:
        for identifier, msa in tandem_repeats.items():
            f.write(">{0}\n".format(identifier))
            f.write("\n".join(msa) + "\n\n")
