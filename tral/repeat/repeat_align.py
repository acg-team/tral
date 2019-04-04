# (C) 2015 Elke Schaper

from tral import configuration
"""
    :synopsis: Alignment of tandem repeat units.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import os
import subprocess
import tempfile

log = logging.getLogger(__name__)


CONFIG_GENERAL = configuration.Configuration.instance().config
REPEAT_CONFIG = CONFIG_GENERAL["repeat"]

''' Some functions might overlap with repeat.gene_tree.align.'''


def realign_repeat(my_msa, aligner='mafft', sequence_type='AA', begin=None, rate_distribution='Constant'):

    # TODO: should we give the possibility to set gamma distribution easily for one self?
    # TODO: rate_distribution=gamma or rate_distribution=gamma(n=6,alpha=0.5) or gamma=TRUE/FALSE?

    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    log.debug("evolvedTR: Created temp directory: %s", working_dir)

    # Save my_TR to temp directory:
    msa_file = os.path.join(working_dir, 'initial_msa.faa')
    with open(msa_file, 'w') as msa_filehandle:
        for i, iMSA in enumerate(my_msa):
            msa_filehandle.write('>t{0}\n{1}\n'.format(i, iMSA))

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
        
    elif aligner == "proPIP":
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
            print("Sequence type is not known.")

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
            print("For two units which have to be aligned an arbritary tree will be given.")            
            tree_string = "(t0:0.1,t1:0.1);"
            init_tree = "user"
            with open(tree_initial, 'w') as treefile:
                treefile.write(tree_string)

        # For three units an arbritary initial tree is used for the alignment
        if len([1 for line in open(msa_file) if line.startswith(">")]) == 3:
            print("For three units which have to be aligned an arbritary tree will be given.")            
            tree_string = "((t0:0.1,t2:0.1):0.1,t1:0.1);"
            init_tree = "user"
            with open(tree_initial, 'w') as treefile:
                treefile.write(tree_string)
        # For more than tree units an initial tree will be estimated from the previous alignment 

        # create parameter file to create tree with castor
        parameters_tree =  ["analysis_name=tree_optimization",
                            "alphabet={}".format(alphabet),
                            "alignment=false",
                            "input.sequence.file={}".format(msa_file),
                            "input.sequence.sites_to_use=all",
                            "init.tree={}".format(init_tree),
                            "input.tree.file={}".format(tree_initial),
                            "init.distance.method=bionj",
                            "model=PIP(model={}(initFreqs=observed),initFreqs=observed)".format(substitution_model), 
                            "rate_distribution={}".format(rate_distribution),
                            "optimization=D-BFGS(derivatives=BFGS)",
                            "optimization.max_number_f_eval=500", 
                            "optimization.tolerance=0.001", 
                            "optimization.final=bfgs",
                            "optimization.topology=true",
                            "optimization.topology.algorithm=Mixed(coverage=best-search,starting_nodes=Hillclimbing(n=4),max_cycles=50,tolerance=0.01,brlen_optimisation=BFGS,threads=10)",
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
                                                            "params={}".format(paramsfile_tree)])
            castor_tree_initialization.wait()

            # TODO: Catch this error: Column #33 of the alignment contains only gaps. Please remove it and try again!
            # the given alignment cannot have a column with only gaps
        except FileNotFoundError:
            error_note = (
                "ProPIP could not be reached.\n" +
                "Is Castor installed properly and is the path defined in config.ini in the data directory?\n")
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
                                "model=PIP(model={},lambda=10,mu=0.5)".format(substitution_model), # add to config?
                                # "model=PIP(model={}(initFreqs=observed),initFreqs=observed)".format(substitution_model),  # add to config?
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
                                                "params={}".format(paramsfile_alignment)])
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
                tree_output = f.readlines()
        except FileNotFoundError:
            try:    
                with open(msa_realigned, "r") as f:
                    tree_output = f.readlines()
            except FileNotFoundError:
                error_note = (
                    "ProPIP could not successfully be used for the realignment of:\n" +
                    "\n".join(my_msa))
                logging.error(error_note)
                return
                
        msa = [iLine[:-1] for iLine in tree_output if iLine[0] != '>']
        label = [iLine[:-1] for iLine in tree_output if iLine[0] == '>']
        # use original order of msa
        msa_sorted = [x for _,x in sorted(zip(label,msa))]

        log.debug('\n'.join(msa_sorted))
        try:
            return msa_sorted
        except:
            error_note = (
                "ProPIP could not successfully be used for the realignment of:\n" +
                "\n".join(my_msa))
            logging.error(error_note)
            return None

    else:
        raise ValueError(
            'Currently, the aligner {} is not implemented.'.format(aligner))
