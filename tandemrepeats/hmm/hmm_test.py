import itertools, logging

from . import hmm, sequence, viterbi
from ..repeat import repeat_info, repeat_io

logger = logging.getLogger('root')
logger.setLevel(logging.WARNING)
#logging.basicConfig(level=logging.WARNING)

def create_test_set(n = 3, l = [3,15,30], divergence = [0], iterations = 10, flank_length_multiplier = 5, sequence_type = 'AA'):

    ''' Create a test set for HMMs consisting of
    <iterations> pairs of evolved homologous tandem repeats of repeat unit length <l> and number of repeat units <n>
    and both-sided random flanks of length <flank_length_multiplier>*<l>.'''

    tests = [{'n': n, 'l': iL, 'n_samples': iterations, 'flank_length':  flank_length_multiplier*iL, 'divergence': iD, 'seed': {}, 'detectable': {}} for iL,iD in itertools.product(l,divergence)]
    
    
    for iT in tests:
        
        # Create Tandem Repeats
        if iT['divergence'] > 0:
            tandem_repeats = list(repeat_io.evolved_tandem_repeats(l=iT['l'], n=2*iT['n'], n_samples=iT['n_samples'], mutationRate=iT['divergence'], sequence_type=sequence_type, return_type = 'list'))
        else:
            tandem_repeats = []
            while len(tandem_repeats) < iT['n_samples']:
                new_repeats = list(repeat_io.random_sequence(n_samples = iT['n_samples']-len(tandem_repeats), sequence_length = iT['l'], return_type = 'list')) 
                tandem_repeats.extend([[iTR]*(2*iT['n']) for iTR in new_repeats if not 'X' in iTR])
            
        
        # Split into seed and detectable tandem repeats
        iT['seed']['tandem_repeats'] = [repeat_info.Repeat(begin = 0, msa = iMSA[:iT['n']]) for iMSA in tandem_repeats]
        iT['detectable']['tandem_repeats'] = [repeat_info.Repeat(begin = 0, msa = iMSA[iT['n']:]) for iMSA in tandem_repeats]
        
        # Create Random Flanks
        random_flanks = []
        while len(random_flanks) < 2*iT['n_samples']:
            new_flanks = list(repeat_io.random_sequence(n_samples = 2*iT['n_samples']-len(random_flanks), sequence_length = 2*iT['flank_length'], return_type = 'list')) 
            random_flanks.extend([iF for iF in new_flanks if not 'X' in iF])
        #random_flanks = list(repeat_io.random_sequence(n_samples = 2*iT['n_samples'], sequence_length = 2*iT['flank_length'], return_type = 'list'))
        
        # Create complete test (seed and detectable) sequences from the first and second half, respectively, of the random flanks set
        iT['seed']['sequences'] = [iF[:iT['flank_length']] + "".join(iTR.msa) + iF[iT['flank_length']:]   for iF,iTR in zip(random_flanks[:iT['n_samples']], iT['seed']['tandem_repeats'] )] 
        iT['detectable']['sequences'] = [iF[:iT['flank_length']] + "".join(iTR.msa) + iF[iT['flank_length']:]   for iF,iTR in zip(random_flanks[iT['n_samples']:], iT['detectable']['tandem_repeats'] )]
            
            
    return tests
    

def run_hmms_on_test_set(tests, divergence, parameter_ID):

    # Create a sequence.Sequence object.
    my_sequence = sequence.Sequence()    
    
    ## Define some kind of cut-off for which repeats to use and which not:
    #cut_off_constraint = {'alpha' : 0.01, 'classifier' : 'phylo_gap01_ignore_trailing_gaps_and_coherent_deletions'}
    
    ## Define classifiers of interest:
    #classifiers = ['phylo', 'phylo_gap01', 'phylo_gap001']

    for iT in tests:
        
        iT[parameter_ID] = {'tandem_repeats': [], 'left_greediness': [], 'right_greediness': []}
        
        for iS in range(iT['n_samples']):
            
            tandem_repeat = iT['seed']['tandem_repeats'][iS]
            
            iHMM = hmm.HMM(tandem_repeat, divergence = divergence) 
    
            # Copy the amino acid sequence for each leaf into the sequence.Sequence object <my_sequence>.
            my_sequence.sequence = iT['detectable']['sequences'][iS]
            
            # Create a viterbi.Viterbi object from the search sequence <my_sequence> and the HMM <my_HMM>.
            my_viterbi = viterbi.Viterbi(iHMM, my_sequence.sequence)
            
            # Run the Viterbi algorithm.
            my_most_likely_path = my_viterbi.viterbi()
            
            # Convert the Viterbi path into a tandem repeat. Calculate the shift towards the original tandem repeat (my_tandem_repeat.shift).
            my_tandem_repeat = viterbi.hmm_path_to_tandem_repeat(my_sequence.sequence, most_likely_path = my_most_likely_path, lD=tandem_repeat.lD)
            
            iT[parameter_ID]['tandem_repeats'].append(my_tandem_repeat)
            
            # If a tandem repeat was found...
            if my_tandem_repeat == None:
                iT[parameter_ID]['left_greediness'].append(None)
                iT[parameter_ID]['right_greediness'].append(None)
                print('no!')
            else:
                iT[parameter_ID]['left_greediness'].append(iT['flank_length'] - my_tandem_repeat.begin)
                iT[parameter_ID]['right_greediness'].append(my_tandem_repeat.begin + my_tandem_repeat.sequence_length - iT['flank_length'] - iT['n']*iT['l'])
            
                
    return tests
    
def evaluate_hmm(tests, divergence = [0.01,0.1,0.5,1]):

    for iD in divergence:
        tests = run_hmms_on_test_set(tests,iD,parameter_ID=iD)        
    return tests