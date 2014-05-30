# (C) 2013 Elke Schaper

import logging, os, subprocess, tempfile
from Bio import AlignIO

logger = logging.getLogger('root')

from tandemrepeats.repeat import repeat_info
from tandemrepeats.paths import *

SCORESLIST = ['phylo_gap01']

''' Some functions might overlap with repeat.gene_tree.align.'''

def realign_repeat(my_msa, aligner = 'mafft', sequence_type = 'AA', begin = None):

    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    logger.debug("evolvedTR: Created temp directory: %s", working_dir)

    # Save my_TR to temp directory:
    msa_file = os.path.join(working_dir, 'msa_temp.faa')
    with open(msa_file, 'w') as msa_filehandle:
        for i,iMSA in enumerate(my_msa):
            msa_filehandle.write('> {0}\n{1}\n'.format(i, iMSA))

    if aligner == 'mafft':
        # Run Mafft
        # See http://mafft.cbrc.jp/alignment/software/manual/manual.html for choice of options.
        # The mafft result is in stdout. Check: Do you need to capture or redirect the stderr?
        p = subprocess.Popen(["ginsi", "--anysymbol", "--quiet", msa_file], stdout=subprocess.PIPE)
        mafft_output = [line.decode('utf8').rstrip() for line in p.stdout]
        msa = []
        for iLine in mafft_output:
            if iLine[0] == '>':
                msa.append('')
            else:
                msa[-1] += iLine
        msa = [i for i in msa if i != '']
        logger.debug('\n'.join(msa))
        p.wait()
        try:
            return msa
            #return repeat_info.Repeat(begin=begin, msa = msa, sequence_type = sequence_type, calc_score = True, calc_pValue = True, scoreslist=SCORESLIST)
        except:
            print("\n".join(my_msa))
            print('*'*3)
            print("\n".join(msa))
            return None

    elif aligner == 'prograph':
        print('Haha, not yet.')


########################################### MAIN #########################################

if __name__=="__main__":

    msa = ["MGKGYL---------------------------------------ALCSYNCKEA-INILSHLPSHHYN","TG--------------------------------------------------------------WVLCQ","IGRAYF---------------------------------------ELSEYMQAER-IFSEVRRIENYRV","EGMEIYSTTLWHLQK------------------------------DVALSVLSKDLTDMDKNSPEAWCA","AGNCFS---------LQREH-------------------------DIAIKFFQRA-IQVDPNYAYAYTL","LGHEFV--------------LTEEL--------------------DKALACFRNA-IRVNPRHYNAWYG","LGMIYY-------------------KQEKF---------------SLAEMHFQKA-LDINPQSSVLLCH","IGVVQH------------------------ALKKS----------EKALDTLNKA-IVIDPKNPLCKFH","RASVLF-----------------------------ANEKY-----KSALQELEEL-KQIVPKESLVYFL","IGKVYK----------------------------------KLGQTHLALMNFSWA-MDLDPKGAN----"]
    realign_repeat(msa = msa)


