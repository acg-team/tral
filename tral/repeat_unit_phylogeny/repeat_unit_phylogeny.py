# (C) 2015 Elke Schaper

"""
    :synopsis: A phylogeny of tandem repeats units

    .. moduleauthor:: Elke Schaper <elke.schaper@isb.sib.ch>
"""

import itertools

import datetime
import logging
import numpy as np
import os
import pickle
import shutil
import subprocess
import tempfile

from ete3 import Tree

from tral.hmm import hmm_viterbi
from tral.repeat.repeat import Repeat
from tral.repeat import repeat_align
from tral import configuration

LOG = logging.getLogger(__name__)

CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat_unit_phylogeny"]

##################### Repeat unit phylogeny class #############################


class RepeatUnitPhylogeny:

    """ Repeat unit phylogenies display the evolution of a tandem repeat, or
        several homolog tandem repeat in terms of unit duplications and
        deletions.

    Note:
        The repeat unit phylogenies implemented here were introduced in:

        Schaper, E., Gascuel, O. & Anisimova, M. Deep conservation of human
        protein tandem repeats within the eukaryotes.
        Molecular Biology and Evolution (2014).



    Attributes:
        id (str): The id describes the origin of the RepeatUnitPhylogeny.
        phylogeny (ete2 phylogeny): The phylogeny

    """

    def __str__(self):
        """
            Create string for RepeatUnitPhylogeny instance.
        """
        return str(self.phylogeny)


    def __init__(self, phylogeny, id = None):
        """ HMM class __init_ module.

        __init__ takes HMM parameters (including the alphabet, emission
        probabilities and transition probabilities) as input, and assigns
        them to class attributes.

        Args:
            hmmer_probabilities (dict): A dictionary with HMM parameters.

        ..todo:: phylogeny may currently both be a string and a file. Implement
                appropriate test.
        """
        self.id = id
        LOG.debug("The phylogeny's ID is {}".format(self.id))

        self.phylogeny = Tree(phylogeny)


    def create(input_type, file=None, repeats=None, **kwargs):
        """ Creates a `RepeatUnitPhylogeny` instance from several input types.

        A `RepeatUnitPhylogeny` instance is created from one of the following
        input types:

        * Single or multiple ``Repeat`` instances
        * files (e.g. newick format)

        Args:
            input_type (str): The input format. E.g. "repeat", "file"
            file (str): Path to the file containing the phylogeny.
            repeats (lists of Repeat): A list of Repeat objects with an MSA
                that can be transformed into an HMM.

        Returns:
            `RepeatUnitPhylogeny`: An initialized instance of the
                `RepeatUnitPhylogeny` class.

        Raises:
            Exception: if the phylogeny file does not exist.
            Exception: if the repeat value provided is not an instance of the
                Repeat class.
            Exception: if no parameters are provided to create the
                RepeatUnitPhylogeny.

        """

        if input_type == 'file':
            # Read phylogeny from file.
            phylogeny = file
        elif input_type == 'repeats':
            if not (isinstance(repeats, Repeat) or isinstance(repeats, list)):
                raise Exception('repeats is neither a valid instance of list '
                                'nor of the Repeat class.')
            elif isinstance(repeats, list) and not all(isinstance(iR, Repeat) for iR in repeats):
                raise Exception('repeats is a list, but the list members are '
                                'not all instances of the Repeat class.')
            if 'PhyML' in kwargs:
                phylogeny = RepeatUnitPhylogeny.create_from_repeat(repeats, **kwargs['PhyML'])
            else:
                phylogeny = RepeatUnitPhylogeny.create_from_repeat(repeats)
            LOG.debug(phylogeny)
        else:
            raise Exception("Unknown input type: {}.".format(input_type))

        return RepeatUnitPhylogeny(phylogeny)


    def create_from_repeat(repeats, result_file=None, method='PhyML', sequence_type='aa'):
        """ Get HMM parameters (including the alphabet, emission probabilities
        and transition probabilities) from a tandem repeat.

        An HMM is created from ``tandem_repeat`` using :command:`hmmbuild`
        and is saven in the HMMER3 file format.
        Next, HMM parameters are retrieved from the HMMER3 file and returned.

        Args:
            tandem_repeat (TR):  A Repeat class instance of the TR to be
                transformed into an HMM.
            tmp_file (str): Path to where a copy of the created phylogeny file
                will be stored. If None, no copies are saved.


        Returns:
            file (str): Path to the phylogeny file.

        """

        # Create temporary working directory.
        tmp_dir = tempfile.mkdtemp()
        LOG.debug("Reconstruct_tree: Created temp directory: %s", tmp_dir)

        # Save tandem_repeats as Phylip file in temporary working directory.
        alignment_file = os.path.join(tmp_dir, 'input_alignment.phylip')
        homologs = {} # homologs e.g. {'ENSP00012': "ABC-E", 'ENSP00013': "ABCDE"}
        for iRepeat in repeats:
            homologs.update({"{}_{}".format(iRepeat.name, i): iUnit for i,iUnit in enumerate(iRepeat.msa)})
        write_msa_phylip(alignment_file, homologs)

        if method == 'PhyML':
            # Run PhyML
            # See http://www.atgc-montpellier.fr/download/papers/phyml_manual_2009.pdf for choice of options.
            # The PhyML result is in stdout?

            #        # Define PhyML sequence type flag.
            #         if tandem_repeat.sequence_type == "AA":
            #             sequence_type_flag = "--amino"
            #         elif tandem_repeat.sequence_type == "RNA":
            #             sequence_type_flag = "--rna"
            #         else:
            #             sequence_type_flag = "--dna"

            # Anything useful in here? :
            #        p = subprocess.Popen([CONFIG["hmmbuild"], sequence_type_flag, tmp_id + ".hmm",
            #                              tmp_id + ".sto"],
            #                             stdout=subprocess.PIPE, stderr=None, cwd=tmp_dir)

            #phyML_path = os.path.join(EXECROOT, 'PhyML_3.0_linux64')
            tree_file = os.path.join(tmp_dir,'input_alignment.phylip_phyml_tree.txt')
            p = subprocess.Popen([CONFIG["PHYML"], "-i", alignment_file, "-d",
                                 sequence_type, "-q", "-s NNI", "-b 0", "--quiet"],
                                stdout=open(os.devnull, 'wb'), close_fds=True)
            p.wait()

            if result_file:
                shutil.copyfile(tree_file, result_file)

            #from Bio import Phylo
            #tree = Phylo.read(tree_file, "newick")

            with open(tree_file,'r') as th:
                tree = th.readline().rstrip()
            LOG.debug(tree)

            shutil.rmtree(tmp_dir)

            return tree


    def write(self, file, output_format, *args):
        """ Write ``RepeatUnitPhylogeny`` to file.

        Write ``RepeatUnitPhylogeny`` to file. Currently, only pickle is
        implemented.

        Args:
            file (str): Path to input file
            output_format (str):  Either "pickle" or "newick".

        .. todo:: Write checks for ``format`` and ``file``.
        .. todo:: Implement further formats.

        """

        if output_format == 'pickle':
            with open(file, 'wb') as fh:
                pickle.dump(self, fh)
        else:
            raise Exception('output_format is unknown: {}'.format(output_format))


def pairwise_optimal_repeat_unit_phylogenies(d_sequences, d_viterbi_path, lD):
    """ Calculate pairwise repeat unit phylogenies for a set of homolog repeats,
        for which only the Viterbi paths (from a cpHMM detection) and the
        sequences are known.

        * Retrieve optimal alignments for each pair
        * Realign pairwise MSAs
        * Build phylogeny of all realigned pairwise repeat units
        * Perhaps: Create .pdf of phylogeny with scriptree.

        Args:
            d_sequences (dict of str): Sequences as str and their identifiers.
            d_viterbi_path (dict of list):  Viterbi paths as lists of str and their identifiers.
            lD (int): Length of the cpHMM used to detect the repeats.

        """

    l_id = list(d_sequences.keys())

    results = {(iV,jV): None for iV,jV in itertools.combinations(l_id, 2)}
    for v in results.keys():

        # Retrieve optimal alignments for each pair
        sequences = [d_sequences[v[0]], d_sequences[v[1]]]
        viterbi_paths = [d_viterbi_path[v[0]], d_viterbi_path[v[1]]]
        l_msa = hmm_viterbi.hmm_path_to_maximal_complete_tandem_repeat_units(sequences, viterbi_paths, lD)

        # Realign pairwise MSAs.
        unaligned_msa = [iUnit for iMSA in l_msa for iUnit in iMSA]
        aligned_msa = repeat_align.realign_repeat(unaligned_msa)

        # Build repeat unit phylogeny. First, make two repeats
        repeat1 = Repeat(msa=aligned_msa[:len(unaligned_msa[0])], name=v[0])
        repeat2 = Repeat(msa=aligned_msa[len(unaligned_msa[0]):], name=v[1])

        # Next, create the phylogeny.
        phylo = RepeatUnitPhylogeny.create(input_type = "repeats", repeats = [repeat1, repeat2])

        # Calculate conservation
        LOG.warning("Not implemented yet.")

        # Make drawings
        LOG.warning("Not implemented yet.")
        #if draw_phylogeny:
        #    # Create scripttree files to later create script trees.
        #    tree_draw_file = os.path.join(result_dir, "{0}".format(key))
        #    #tree_io.create_scriptree_files(result_file_stump = os.path.join(result_dir, "{}_{}_{}".format(my_TR_ID,ensembl_ID,iO)), tree = my_TR_data['homologous_pairs'][ensembl_ID][iO], annotations = None)
        #    tree_io.create_scriptree_files(tree_draw_file, tree_nhx = newick_tree, annotations = annotations)

        results[v] = phylo

    return results


def write_msa_phylip(file, homologs):
    ''' Write multiple homologs in sequential PHYLIP format to `file`.

        Write multiple homologs in sequential PHYLIP format to `file`. Example
        of the PHYLIP format:
        2 5
        ENSP00012 ABC-E
        ENSP00013 ABCDE

        Args:
            file (str): Path to input file
            homologs (dict of str):  dict of aligned homologs and identifiers,
                e.g. {'ENSP00012': "ABC-E", 'ENSP00013': "ABCDE"}


        ..todo:: Refactor method to `repeat`
    '''

    if len(homologs) == 0:
        raise ValueError('Dict <homologs> is empty.')
        return None

    with open(file, 'w', newline = '\n') as f:

        f.write("{0} {1}\n".format(len(homologs), len(next(iter(homologs.values()))) ))
        for identifier,aligned_sequence in homologs.items():
            f.write("{0} {1}\n".format(identifier, aligned_sequence))

    return True
