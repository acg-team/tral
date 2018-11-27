# (C) 2015 Elke Schaper

"""
    :synopsis: Alignment of tandem repeat units.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import os
import subprocess
import tempfile
from Bio import AlignIO

log = logging.getLogger(__name__)

from tral.repeat import repeat
from tral import configuration

CONFIG_GENERAL = configuration.Configuration.instance().config
REPEAT_CONFIG = CONFIG_GENERAL["repeat"]

''' Some functions might overlap with repeat.gene_tree.align.'''


def realign_repeat(my_msa, aligner='mafft', sequence_type='AA', begin=None):

    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    log.debug("evolvedTR: Created temp directory: %s", working_dir)

    # Save my_TR to temp directory:
    msa_file = os.path.join(working_dir, 'msa_temp.faa')
    with open(msa_file, 'w') as msa_filehandle:
        for i, iMSA in enumerate(my_msa):
            msa_filehandle.write('> {0}\n{1}\n'.format(i, iMSA))

    if aligner == 'mafft':
        # Run Mafft
        # See http://mafft.cbrc.jp/alignment/software/manual/manual.html for choice of options.
        # The mafft result is in stdout. Check: Do you need to capture or
        # redirect the stderr?
        p = subprocess.Popen([REPEAT_CONFIG['ginsi'],
                              "--anysymbol",
                              "--quiet",
                              msa_file],
                             stdout=subprocess.PIPE)
        mafft_output = [line.decode('utf8').rstrip() for line in p.stdout]
        msa = [iLine for iLine in mafft_output if iLine[0] != '>']
        log.debug('\n'.join(msa))
        p.wait()
        try:
            return msa
        except:
            error_note = (
                "Mafft could not successfully run the realignment for: "
                "\n".join(my_msa))
            logging.error(error_note)
            return None

    else:
        raise ValueError(
            'Currently, the aligner {} is not implemented.'.format(aligner))
