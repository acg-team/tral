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


def realign_repeat(my_msa, aligner='mafft', sequence_type='AA', begin=None):

    # Create temporary working directory
    working_dir = tempfile.mkdtemp()
    log.debug("evolvedTR: Created temp directory: %s", working_dir)

    # Save my_TR to temp directory:
    msa_file = os.path.join(working_dir, 'msa_temp.faa')
    with open(msa_file, 'w') as msa_filehandle:
        for i, iMSA in enumerate(my_msa):
            msa_filehandle.write('>{0}\n{1}\n'.format(i, iMSA))

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

        # TODO: change substitution models
        # TODO: put alignment into std output instead of file
        # TODO: include gamma option
        # TODO: Use correct naming for Castor/proPIP. proPIP should not be referred to as castor!!

        # log messages of castor to stderr instead of logfiles
        # os.environ["GLOG_logtostderr"] = "1"

        if sequence_type == "AA":
            alphabet="Protein"
        elif sequence_type == "DNA":
            alphabet="DNA"
        else:
            "Sequence type is not known."

        tree = os.path.join(working_dir,"tree.nwk")
        msa_realigned = os.path.join(working_dir,"msa_realigned.faa")

        ####################################
        # Create an initial tree
        ####################################

        # Castor cannot create trees for less than four sequences

        # Aligning two units is currently not supported
        if len([1 for line in open(msa_file) if line.startswith(">")]) == 2:
            raise Exception('Sorry, a tandem repeat with only two sequences cannot be realigned with proPIP.')

        # For three units an arbritary initial tree is used for the alignment
        elif len([1 for line in open(msa_file) if line.startswith(">")]) == 3:
                print("For three units which have to be aligned an arbritary tree will be given.")            
                tree_string = "((0:0.1,2:0.1):0.1,1:0.1);"
                with open(tree, 'w') as treefile:
                    treefile.write(tree_string)

        # For more than tree units an initial tree will be estimated from the previous alignment 
        else:
            try:
                castor_tree_initialization = subprocess.Popen([REPEAT_CONFIG['Castor'],
                                    "analysis_name=prova",
                                    "model_description=JC69+PIP",  # add to config?
                                    "input_folder={}".format(working_dir),
                                    "output_folder={}".format(working_dir),
                                    "alphabet={}".format(alphabet),
                                    "alignment=false",
                                    "input.sequence.file={}".format(msa_file),
                                    "input.sequence.sites_to_use=all",
                                    "init.tree=distance",
                                    "init.distance.method=bionj",
                                    "model=PIP(model=JC69(initFreqs=observed),initFreqs=observed)",  # add to config?
                                    "rate_distribution=Constant",
                                    "optimization=D-BFGS(derivatives=BFGS)",
                                    "optimization.max_number_f_eval=5000",  # add to config?
                                    "optimization.tolerance=0.001",  # add to config?
                                    "optimization.final=bfgs",
                                    "optimization.topology=true",
                                    "optimization.topology.algorithm=Mixed(coverage=best-search,starting_nodes=Hillclimbing(n=8),max_cycles=100,tolerance=0.001,brlen_optimisation=BFGS,threads=10)",  # add to config?
                                    "output.estimates.file={}".format(working_dir),
                                    "output.tree.file={}".format(tree),
                                    "output.estimates.format=json"
                                    "support=none"])
                castor_tree_initialization.wait() # needed that the next analysis can be executed correctly! 
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
            print("A problem occurred while trying to reach previous alignment in file.") # TODO: Improve this error message!

        try:
            proPIP_alignment = subprocess.Popen([REPEAT_CONFIG['Castor'],
                                "analysis_name=aligner",
                                "model_description=GTR+PIP",  # add to config?
                                "input_folder={}".format(working_dir),
                                "output_folder={}".format(working_dir),
                                "alphabet={}".format(alphabet),
                                "alignment=true",
                                "alignment.version=ram",
                                "input.sequence.file={}".format(unaligned_sequences),
                                "input.sequence.sites_to_use=all",
                                "init.tree=user",
                                "init.distance.method=bionj",
                                "input.tree.file={}".format(tree),
                                "model=PIP(model=JC69,lambda=0.2,mu=0.1)", # add to config?
                                "rate_distribution=Constant",
                                "optimization=None",
                                "output.msa.file={}".format(msa_realigned),
                                "output.estimates.file={}".format(working_dir),
                                "output.estimates.format=json"
                                "support=none"])
            proPIP_alignment.wait()
        except FileNotFoundError:
            error_note = (
                "ProPIP could not be reached.\n" +
                "Is Castor installed properly and is the path defined in config.ini in the data directory?\n")
            logging.error(error_note)
            return

        try:
            # The created alignment file has "initial" included into the name because in future the tool should be able to realign
            with open(os.path.join(working_dir,"msa_realigned.initial.faa"), "r") as f:
                castor_output = f.readlines()
        except FileNotFoundError:
            try:    
                with open(msa_realigned, "r") as f:
                    castor_output = f.readlines()
            except FileNotFoundError:
                error_note = (
                    "Castor could not successfully run the realignment for:\n" +
                    "\n".join(my_msa))
                logging.error(error_note)
                return
                
        msa = [iLine[:-1] for iLine in castor_output if iLine[0] != '>']
        label = [iLine[:-1] for iLine in castor_output if iLine[0] == '>']
        # use original order of msa
        msa_sorted = [x for _,x in sorted(zip(label,msa))]

        log.debug('\n'.join(msa_sorted))
        try:
            return msa_sorted
        except:
            error_note = (
                "Castor could not successfully run the realignment for:\n" +
                "\n".join(my_msa))
            logging.error(error_note)
            return None

    else:
        raise ValueError(
            'Currently, the aligner {} is not implemented.'.format(aligner))
