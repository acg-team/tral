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
        # print("proPIP is not integrated yet.")
        print("Start doing something with Castor.")

        if sequence_type == "AA":
            alphabet="Protein"
        elif sequence_type == "DNA":
            alphabet="DNA"
        else:
            "Sequence type is not known."

        tree = os.path.join(working_dir,"tree.nwk")
        msa_realigned = os.path.join(working_dir,"msa_realigned.faa")

        # create file with unaligned sequences
        unaligned_sequences = os.path.join(working_dir,"sequences.faa")                    
        with open(msa_file, 'r') as infile, open(unaligned_sequences, 'w') as outfile:
            temp = infile.read().replace("-", "")
            outfile.write(temp)

        # TODO:
        # - replace Castor with config-castor
        # - replace files in parameter file
        # - prevent printing of output files in current directory!

        ###################
        # add castor to the config file

        ###################
        # call castor to infer a tree from the previous alignment
        

        castor_tree_initialization = subprocess.Popen(["Castor",
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

        # call castor to align TR units with inferred tree 

        # TODO: model dependent on DNA/Protein

        castor_alignment = subprocess.Popen(["Castor",
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
        castor_alignment.wait()

        # TODO: put alignment into std output instead of file
        # TODO: change name of alignment! this initial makes it strange....


        with open(os.path.join(working_dir,"msa_realigned.initial.faa"), "r") as f:
            castor_output = f.readlines()

        msa = [iLine[:-1] for iLine in castor_output if iLine[0] != '>']
        log.debug('\n'.join(msa))
        try:
            return msa
            # TODO: see if the msa has to be in a correct turn... if yes sort the msa file before creating list
        except:
            error_note = (
                "Castor could not successfully run the realignment for: "
                "\n".join(my_msa))
            logging.error(error_note)
            return None



    else:
        raise ValueError(
            'Currently, the aligner {} is not implemented.'.format(aligner))
