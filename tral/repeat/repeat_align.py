# (C) 2015 Elke Schaper


"""
    :synopsis: Alignment of tandem repeat units.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
    .. moduleauthor:: Paulina Näf
"""

from tral import configuration
import logging
import os
import subprocess
import tempfile

log = logging.getLogger(__name__)

CONFIG_GENERAL = configuration.Configuration.instance().config
REPEAT_CONFIG = CONFIG_GENERAL["repeat"]

def realign_repeat(my_msa, realignment='mafft', sequence_type='AA', rate_distribution=REPEAT_CONFIG['castor_parameter']['rate_distribution'], user_path=None):

    """ Realignment of a repeat MSA using mafft or proPIP

    A good multiple sequence alignment (MSA) is an important description of a
    tandem repeat (TR) and can be crucial to distinguish a true positive from
    a false positive TR. 
    This can be realised using a proper MSA estimation algorithm. Currently,
    Mafft and proPIP can be used to realign a TR.

    For proPIP the rate distribution for indels can be either constant or gamma distributed
        - REPEAT_CONFIG['castor_parameter']['rate_distribution'] == 'constant': no gamma distribution is used for the proPIP algorithm
        - REPEAT_CONFIG['castor_parameter']['rate_distribution'] == 'gamma': gamma distribution (n=6, alpha=0.5)

    Args: 
        my_msa (list of strings): List of sequences (str)
        realignment (str): Either "mafft" or "proPIP"
        sequence_type (str): Either "AA" or "DNA"
        user_path (str): copy alignment files to user defined path
    
    Returns:
        my_msa realigned (list of strings)

    .. todo:: 
        - test if input alignment has column with only gaps
        - decide what happens if branch length of tree is getting zero
    """

    realignment_types = ['mafft', 'proPIP', None]
    if realignment not in realignment_types:
        raise ValueError("Invalid realignment option. Expected one of: {}".format(realignment_types))
    
    if not all(isinstance(s, str) for s in my_msa):
        raise ValueError("Invalid MSA input.")

    if realignment == "proPIP":
        # remove columns that only contains gaps to avoid errors
        my_msa = remove_gaps(my_msa)

        # do try to realign if no gap is present in MSA
        if not any(['-' in s for s in my_msa]):
            print("\nproPIP needs at least one gap in the MSA to realign.\n" +
                "No realignment will be performed.\n")
            return my_msa

    if realignment == None:
        print("\nNo realignment was performed.")
        return my_msa

    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    log.debug("evolvedTR: Created temp directory: %s", working_dir)

    # Save my_TR to temp directory:
    msa_file = os.path.join(working_dir, 'initial_msa.faa')
    with open(msa_file, 'w') as msa_filehandle:
        for i, iMSA in enumerate(my_msa):
            msa_filehandle.write('>{0}\n{1}\n'.format(i+1, iMSA))

    if realignment == 'mafft':
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
        msa = []
        label= []
        for i, line in enumerate(mafft_output):
            if line[0] == '>':
                label.append(int(line[1:]))
            elif mafft_output[i-1][0] == '>':
                msa.append(line)
            else:
                msa[-1] += line
        # use original order of msa
        msa_sorted = [x for _,x in sorted(zip(label,msa))]

        log.debug('\n'.join(msa_sorted))
        p.wait()
        try:
            return msa
        except:
            error_note = (
                "Mafft could not successfully run the realignment for: "
                "\n".join(my_msa))
            logging.error(error_note)
            return None
        
    elif realignment == 'proPIP':
        # Run Castor (with integrated aligner) (https://github.com/acg-team/castor_aligner)

        # log messages of castor to stderr instead of logfiles
        os.environ["GLOG_logtostderr"] = "1"

        if sequence_type == "AA":
            alphabet="Protein"
            substitution_model = "LG08"
        elif sequence_type == "DNA":
            alphabet="DNA"
            substitution_model="HKY85"
        else:
            raise ValueError("Sequence type is not known.")

        if rate_distribution == 'constant':
            rate_distribution = 'Constant'
        elif rate_distribution == 'gamma':
            rate_distribution = 'Gamma(n=6,alpha=0.5)'
        else:
            raise ValueError("Rate distribution parameter not known.")

        tree_initial = os.path.join(working_dir,"tree_initial.nwk")
        tree = os.path.join(working_dir,"tree.nwk")
        msa_realigned = os.path.join(working_dir,"msa_realigned.faa")
        paramsfile_tree = os.path.join(working_dir, 'params_tree.txt')
        paramsfile_alignment = os.path.join(working_dir, 'params_alignment.txt')
        estimates = os.path.join(working_dir, 'estimates.json')

        ####################################
        # Create an initial tree
        ####################################

        init_tree = "distance"
        # Castor cannot create trees for less than four sequences

        if len([1 for line in open(msa_file) if line.startswith(">")]) == 2:
            print("For two units which have to be aligned an arbritary starting tree will be given.")            
            tree_string = "(1:0.1,2:0.1);"
            init_tree = "user"
            with open(tree_initial, 'w') as treefile:
                treefile.write(tree_string)

        # For three units an arbritary initial tree is used for the alignment
        if len([1 for line in open(msa_file) if line.startswith(">")]) == 3:
            print("For three units which have to be aligned an arbritary starting tree will be given.")            
            tree_string = "((1:0.1,3:0.1):0.1,2:0.1);"
            init_tree = "user"
            with open(tree_initial, 'w') as treefile:
                treefile.write(tree_string)
        # For more than tree units an initial tree will be estimated from the previous alignment 

        # create parameter file to create or optimize tree with castor
        parameters_tree =  ["analysis_name=tree_optimization",
                            "alphabet={}".format(alphabet),
                            "alignment=false",
                            "input.sequence.file={}".format(msa_file),
                            "input.sequence.sites_to_use=all",
                            "init.tree={}".format(init_tree),
                            "input.tree.file={}".format(tree_initial),
                            "init.distance.method=bionj",
                            "model=PIP(model={}(initFreqs=observed),initFreqs=observed)".format(substitution_model), 
                            "rate_distribution=Constant", # creating tree alwas with constant rates to estimate lamda and mu from data
                            "optimization=D-BFGS(derivatives=BFGS)",
                            "optimization.max_number_f_eval=500", 
                            "optimization.tolerance=0.001", 
                            "optimization.final=bfgs",
                            "optimization.topology=true",
                            "optimization.topology.algorithm=Swap(coverage=nnr-search,starting_nodes=Hillclimbing(n=4),max_cycles=50,tolerance=0.01,brlen_optimisation=Brent,threads=10)",
                            "output.tree.file={}".format(tree),
                            "output.tree.format=Newick",
                            "output.estimates.format=json",
                            "output.estimates.file={}".format(estimates),
                            "support=none"]
        try:
            with open(paramsfile_tree, 'w') as params:
                for parameter in parameters_tree:
                    params.write(parameter + '\n')
        except:
            print("A problem occurred while trying to write alignment parameters to txt file")
        
        # run castor for tree initialization
        try:
            castor_tree_initialization = subprocess.Popen([REPEAT_CONFIG['Castor'], 
                                                            "params={}".format(paramsfile_tree)],
                                                            stdout=subprocess.PIPE,
                                                            stderr=subprocess.PIPE)
            print("Tree initialization with Castor ....")
            castor_tree_initialization.wait()

        except FileNotFoundError:
            error_note = (
                "ProPIP could not be reached.\n" +
                "Is Castor installed properly and is the path defined correctly in config.ini in the data directory?\n")
            logging.error(error_note)
            raise Exception("Sorry, creating a tree for the alignment went wrong.")

        ###################
        # Use proPIP algorithm to align TR units with inferred tree
        ###################

        # create file with unaligned sequences
        unaligned_sequences = os.path.join(working_dir,"sequences.faa")
        try:                
            with open(msa_file, 'r') as infile, open(unaligned_sequences, 'w') as outfile:
                temp = infile.read().replace("-", "")
                outfile.write(temp)
        except:
            print("A problem occurred reading initial alignment file.")

        # parse estimates file to use calculated indel parameters
        try:
            import json
            with open(estimates) as estimates_file:
                est = json.load(estimates_file)
                mu_estimated = est['Model']['PIP']['mu']
                lambda_estimated = est['Model']['PIP']['lambda']
                indel_parameters = ',lambda={},mu={}'.format(lambda_estimated,mu_estimated)
        except FileNotFoundError as error:
            print('\nNo file with parameter estimates found.')
            print('Probably it went something wrong when trying to infer or optimise the phylogenetic tree.')
            logging.error(error)
            raise

        # TODO: proPIP cannot produce an alignment when a branchlength is zero.
        #       Need to catch this before an error is thrown
        #       Maybe replace zero with a very small number?
        #       Is this scenario even possible after tree calculation with castor?

        # create parameter file for alignment with proPIP
        parameters_alignment = ["analysis_name=aligner",
                                "model_description={}+PIP".format(substitution_model),
                                "alphabet={}".format(alphabet),
                                "alignment=true",
                                "alignment.version=ram",
                                "input.sequence.file={}".format(unaligned_sequences),
                                "input.sequence.sites_to_use=all",
                                "init.tree=user",
                                "init.distance.method=bionj",
                                "input.tree.file={}".format(tree),
                                "model=PIP(model={}{})".format(substitution_model,indel_parameters),
                                "rate_distribution={}".format(rate_distribution),
                                "optimization=None",
                                "output.msa.file={}".format(msa_realigned),
                                "support=none"]
        try:
            with open(paramsfile_alignment, 'w') as params:
                for parameter in parameters_alignment:
                    params.write(parameter + '\n')
        except:
            print("A problem occurred while trying to write alignment parameters to txt file")
        
        # run alignment with propip
        try:
            proPIP_alignment = subprocess.Popen([REPEAT_CONFIG['Castor'], 
                                                "params={}".format(paramsfile_alignment)],
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)
            print("Realignment of MSA with proPIP .... ")
            proPIP_alignment.wait()
        except FileNotFoundError:
            error_note = (
                "ProPIP algorithm could not be reached.\n" +
                "Is Castor (inkl. proPIP aligner) installed properly and is the path defined in config.ini in the data directory?\n")
            logging.error(error_note)
            return

        # read msa and return as sorted list
        try:
            # The created alignment file has "initial" included into the name because in future the tool should be able to realign
            with open(os.path.join(working_dir,"msa_realigned.initial.faa"), "r") as f:
                realigned = f.readlines()
        except FileNotFoundError:
            try:    
                with open(msa_realigned, "r") as f:
                    realigned = f.readlines()
            except FileNotFoundError:
                error_note = (
                    "ProPIP could not successfully be used for the realignment of:\n" +
                    "\n".join(my_msa))
                logging.error(error_note)
                return
        msa = []
        label = []
        for i, line in enumerate(realigned):
            if line[0] == '>':
                label.append(int(line[1:-1]))
            elif realigned[i-1][0] == '>':
                msa.append(line[:-1])
            else:
                msa[-1] += line[:-1]
        # use original order of msa
        msa_sorted = [x for _,x in sorted(zip(label,msa))]
        log.debug('\n'.join(msa_sorted))

        ## copy alignment files to user defined path
        if user_path:
            if not os.path.exists(os.path.dirname(user_path)):
                os.mkdir(os.path.dirname(user_path))
            import shutil
            mafft_path = user_path.split("proPIP")[0] + "mafft"
            with open(user_path, 'w') as out_msa:
                for i in range(len(msa_sorted)):
                    out_msa.write(">{}\n".format(i+1))
                    out_msa.write("{}\n".format(msa_sorted[i]))
            shutil.copy(msa_file, mafft_path)

        try:
            return msa_sorted
        except:
            raise RuntimeError(
                "ProPIP could not successfully be used for the realignment of:\n" +
                "\n".join(my_msa))
    else:
        raise ValueError(
            'Currently, the aligner {} is not implemented.'.format(realignment))

def remove_char(str, n):
    """ helper function to remove a character in a sequence
    """
    first_part = str[:n] 
    last_part = str[n+1:]
    return first_part + last_part 

def remove_gaps(msa):
    """ Remove columns in a MSA which contain only gaps.
    Do not use MSA longer than 3000 characters.

    Returns the MSA as a list without only-gap-containing columns.

    Args: 
        msa (list of strings): List of sequences (str)
    
    Returns:
        msa (list of strings) with removed colums in case they only contain gaps
    """

    if not all(isinstance(s, str) for s in msa):
        raise ValueError("Invalid MSA input.")

    i = 3000 # The MSA should have at maximum 3000 characters

    if not all(len(s) < i for s in msa):
        raise ValueError("Length of sequences in MSA need to be below {}.".format(i))

    while i >= len(msa[0]):
        i -=1

    while i <= len(msa[0])-1 and not i < 0:
        column = [] 
        for seq in msa:
            column.append(seq[i])
        # print(column)
        if all([char == '-' for char in column]):
            col = 0
            for seq in msa:
                seq = remove_char(seq,i)
                msa[col] = seq
                col += 1
        i -= 1
 
    return msa