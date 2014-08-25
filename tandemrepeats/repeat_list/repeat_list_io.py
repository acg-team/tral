# (C) 2014 Elke Schaper


"""

    :synopsis: Input/output for the repeat_list class

    .. moduleauthor:: Elke Schaper <elke@inf.ethz.ch>

"""
import configobj
import logging
import os

from tandemrepeats.repeat_list.repeat_list import Repeat_list
from tandemrepeats.paths import *

log = logging.getLogger(__name__)

pDefaults = os.path.join(CODEROOT, 'tandemrepeats', 'data', 'defaults.ini')
pSpec = os.path.join(CODEROOT, 'tandemrepeats', 'data', 'spec.ini')
config_general = configobj.ConfigObj(pDefaults, configspec = pSpec)
config = config_general["repeat_list"]


def serialize_repeat_list_tsv(tandem_repeats, config = config):

    ''' Serialize a ``repeat_list`` instance as tsv.

        The following information - if available - is added to the .tsv output:

        * position of the tandem repeats within the sequence,
        * statistical significance of the tandem repeats
        * divergence of the tandem repeat units
        * lengths of the tandem repeat units
        * number of tandem repeat units


    Attributes:
        tandem_repeats (repeat_list): A ``Repeat_list`` instance.
        *args: Additional arguments

    Returns:
        str: The csv as a string.

    '''

    config = {i:j for i,j in config.items() if j}

    tsv = ""

    data = ["\t".join(config.keys())]
    for iRepeat in tandem_repeats.repeats:
        d = []
        for iConfig, type in config.items():
            if "pValue" == iConfig:
                if hasattr(iRepeat,"dPValue") and type in iRepeat.dPValue:
                    d.append( iRepeat.dPValue['type'] )
                else:
                    d.append( None )
            elif "score" == iConfig:
                if hasattr(iRepeat,"dScore") and type in iRepeat.dScore:
                    d.append( iRepeat.dScore['type'] )
                else:
                    d.append( None )
            elif "msa_original" == iConfig:
                try:
                    d.append( ",".join(iRepeat.msa_original) )
                except:
                    raise Exception("The attribute msa_original is not available for iRepeat.")
            elif getattr(iRepeat, iConfig):
                d.append( getattr(iRepeat, iConfig) )
            else:
                raise Exception("The attribute {} is not available for tandem_repeats".format(iConfig))
        data.append( "\t".join(str(i) for i in d) )

    return "\n".join(data)



def save_repeat_fasta(tandem_repeats, file):
    ''' save multiple <tandem_repeats> in Fasta format in specified <file>

        At current, only one TR per sequence can be defined, as the identifiers in
        the dict <tandem_repeats> must be unique.

        Parameters: Dict of tandem repeats and identifiers.
            e.g. {'ENSP00012': msa1, 'ENSP00013': msa2}

            >ID
            GHKI
            GHKI
            GH--
    '''

    with open(file, 'w', newline = '\n') as f:
        for identifier,msa in tandem_repeats.items():
            f.write(">{0}\n".format(identifier))
            f.write("\n".join(msa)+"\n\n")