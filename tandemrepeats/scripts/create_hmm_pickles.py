import argparse
import logging
import logging.config
import os
from Bio import SeqIO

from tandemrepeats.paths import *
from tandemrepeats.hmm import hmm
from tandemrepeats.hmm import hmm_io


logging.config.fileConfig(os.path.join(CODEROOT,'tandemrepeats','data','logging.ini'))
log = logging.getLogger('root')



def main():

    pars = read_commandline_arguments()

    for hmmer_probabilities in hmm_io.read(pars["hmm_file"]):
        iID = hmmer_probabilities['id'].split(".")[0]
        iHMM = hmm.HMM(hmmer_probabilities)
        file = os.path.join(pars['output_path'], "{}.pickle".format(iID))
        iHMM.write(file, "pickle")
        #logging.debug("Wrote {} to {}".format(iHMM.id, pars['output_path']))

def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process create hmm pickles options')
    parser.add_argument('hmm_file', metavar='hmm_file', type=str,
                       help='The path to the input sequence file.')
    parser.add_argument('-o','--output_path', type=str,
                       help='The path to the output files, e.g. /path/to/output')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()
