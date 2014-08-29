#!python

import argparse
import json
import sys

from tandemrepeats import configuration
from tandemrepeats.paths import *
from tandemrepeats.sequence import sequence

c = configuration.Configuration.Instance()
config = c.config

def main():

    pars = read_commandline_arguments()
    print(pars)
    sequence_file = pars.path
    sequence_type = pars.sequence_type
    lTRD = pars.detectors
    dSignificanceTest = json.loads(pars.significance_test)

#     seq = sequence.Sequence.read(fasta_file)[0]
#     test_parameters = {"detection": {"lFinders": ["XSTREAM"]}}
#     denovo_list = seq.detect(denovo = True, **test_parameters)
#     sys.stdout.write(denovo_list.write("tsv") + "\n")

def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('path', metavar='seq.faa', type=str,
                       help='The path to the sequence fasta file.')
    parser.add_argument('-seq','--sequence_type', type=str,
                       help='The sequence type: -seq AA or -seq DNA')
    parser.add_argument('-d', '--detectors', nargs='+', type=str,
                        help='The de novo tandem repeat detectors. For example: -d T-REKS XSTREAM')
    parser.add_argument("-test", "--significance_test", type=str, required=False,
                        help='The significance tests and cut-off values used. For example: -test \'{"phylo_gap01":0.01}\'')
    parser.add_argument("-c", "--cluster", type=str, required=False,
                        help='The overlap definition by which tandem repeats are clustered. For example: -c common_ancestry')
    parser.add_argument("-cf", "--cluster_filter", type=str, required=False,
                        help='The filter by which best representatives are chosen. For example:')


    return parser.parse_args()


if __name__ == "__main__":
    main()