import argparse
import logging
import logging.config
import os
from Bio import SeqIO

from tandemrepeats.paths import *

logging.config.fileConfig(os.path.join(CODEROOT,'tandemrepeats','data','logging.ini'))
log = logging.getLogger('root')


# Code derived from:
#http://biopython.org/wiki/Split_large_file

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = next(iterator)
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch


def main():

    pars = read_commandline_arguments()

    record_iter = SeqIO.parse(open(pars['path']),pars['format'])
    for i, batch in enumerate(batch_iterator(record_iter, pars['number'])) :
        filename = "{}{}.{}".format(pars['output'], str(i+1), pars['format'])
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, pars['format'])
        handle.close()
        logging.debug("Wrote {} records to {}".format(count, filename))
        print("Wrote {} records to {}".format(count, filename))

def read_commandline_arguments():

    parser = argparse.ArgumentParser(description='Process tandem repeat detection options')
    parser.add_argument('path', metavar='sequence_file', type=str,
                       help='The path to the input sequence file.')
    parser.add_argument('-o','--output', type=str,
                       help='The path to the output files, e.g. /path/to/output/name_')
    parser.add_argument('-f', '--format', type=str,
                        help='The file format')
    parser.add_argument("-n", "--number", type=int,
                        help='The number of sequences per file in the output.')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value != None}
    return pars


if __name__ == "__main__":
    main()
