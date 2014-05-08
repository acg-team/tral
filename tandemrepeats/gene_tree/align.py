# (C) 2013 Elke Schaper

import logging, os, shutil, subprocess, tempfile
from Bio import AlignIO

logger = logging.getLogger('root')
from ..repeat import repeat_io
from repeat.paths import *

def read_alignment(file, begin, sequence_type = 'AA'):

    ''' so far, function is not needed. '''

    try:
        alignment = AlignIO.read(file,'fasta')
    except:
        log.error('Alignment output file was not found: {0}'.format(alignment))
        return False
        
    msa = [str(i.seq) for i in alignment]
    # If begin is not known, use repeat_info.repeat_in_sequence function
    return repeat_info.Repeat(begin=0, msa = msa, sequence_type = sequence_type)     
    
def align(sequences, aligner = 'MAFFT', sequence_type = 'AA', tandem_repeat_annotations = None):
    
    '''
    Multiple sequence align <sequences>. If <tandem_repeat_annotations> are given, use <prograph>
    and provide them as added information in t-reks format.
    
    Parameters: Dict of sequences and identifiers
    e.g. {'ENSP00012': 'ABCABAC', 'ENSP00013': 'BCABCA'}   '''
    
    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    logger.debug("align: Created temp directory: %s", working_dir)
    
    if tandem_repeat_annotations:
        aligner = 'ProGraphMSA'
    
    # Save sequences to temp directory:
    sequences_file = os.path.join(working_dir, 'input_sequences.faa')

    ids = []
    with open(sequences_file, 'w') as f:
        for id,sequence in sequences.items():
            f.write('>{0}\n{1}\n'.format(id, sequence))
            ids.append(id)

    if tandem_repeat_annotations:
        tandem_repeat_file = os.path.join(working_dir, 'tandem_repeats.treks')
        tandem_repeat_annotations = {id: iTR for id, iTR in tandem_repeat_annotations.items() if id in ids}     
        # Save homologous tandem repeats to <result_file_treks>
        repeat_io.save_repeat_treks(tandem_repeat_annotations, tandem_repeat_file)
     
    if aligner == 'MAFFT':
        # Run Mafft
        # See http://mafft.cbrc.jp/alignment/software/manual/manual.html for choice of options.
        # The mafft result is in stdout. Check: Do you need to capture or redirect the stderr?
        p = subprocess.Popen(["ginsi", "--anysymbol", "--quiet", sequences_file], stdout=subprocess.PIPE)
        mafft_output = [line.decode('utf8').rstrip() for line in p.stdout]
        msa = {}
        for iLine in mafft_output:
            if iLine[0] == '>':
                id = iLine[1:]
                msa[id] = ''
            else:
                msa[id] += iLine
        logger.debug('\n'.join(msa.values()))
        p.wait()
        return msa

       
    elif aligner == 'ProGraphMSA':                
        # Run ProGraphMSA # Add "--no_force_align" flag if you want Prograph to ignore Methionines in the beginning.
        # Downloaded from http://sourceforge.net/projects/prographmsa/ 
        # The ProGraphMSA result is in stdout. Check: Do you need to capture or redirect the stderr?
        p = subprocess.Popen(["./ProGraphMSA+TR.sh", "--nwdist", "-R","--no_force_align", "-i 0", "--read_repeats", tandem_repeat_file, sequences_file],
                        stdout=subprocess.PIPE, stderr=None, cwd='/cluster/home/infk/eschaper/bin/ProGraphMSA')
        output = [line.decode('utf8').rstrip() for line in p.stdout]
        msa = {}
        for iLine in output:
            if iLine[0] == '>':
                id = iLine[1:]
                msa[id] = ''
            else:
                msa[id] += iLine
                
        if msa == {}:
            shutil.rmtree(os.path.join(ROOT, 'spielwiese'))
            shutil.copytree(working_dir, os.path.join(ROOT, 'spielwiese'))
            raise Exception("ProGraphMSA did not return a MSA. See input files in {0}.".format(os.path.join(ROOT, 'spielwiese')))
                
        #msa = {i for i in msa if i != '']
        logger.debug('\n'.join(msa.values()))
        p.wait()
        return msa
        
def save_msa_fasta(msa, file, tandem_repeat_position = None):

    ''' 
        Parameters <msa>: 
        E.g. {'identifier_1': '----ABCDEF', 'identifier_2': 'GHIJABC---'}
    '''
    
    with open(file,'w') as fh:
        if tandem_repeat_position:
            fh.write("#{0} {1}\n".format(tandem_repeat_position['begin'], tandem_repeat_position['end']))
        for id,seq in msa.items():
            fh.write(">{0}\n{1}\n".format(id, seq))  
            
            
def save_msa_phylip(homologs, file):
    ''' save multiple <homologs> in sequential PHYLIP format in specified <file> 

        Parameters: Dict of aligned homologous sequences and identifiers. 
            e.g. {'ENSP00012': "ABC-E", 'ENSP00013': "ABCDE"}
            
            
        2 5
        ENSP00012 ABC-E
        ENSP00013 ABCDE

    '''
    
    if len(homologs) == 0:
        raise ValueError('Dict <homologs> is empty.')
        return None
    
    with open(file, 'w', newline = '\n') as f:
        
        f.write("{0} {1}\n".format(len(homologs), len(next(iter(homologs.values()))) ))
        for identifier,aligned_sequence in homologs.items():
            f.write("{0} {1}\n".format(identifier, aligned_sequence))
    
    return True 
