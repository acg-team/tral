# (C) 2015 Elke Schaper

"""
    :synopsis: Input/output for tandem repeats.

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import Bio.Seq
import random
import logging
import itertools
import os
import re
import shutil
import tempfile
import subprocess
import numpy as np

from tral.paths import DATA_DIR, EXEC_DIR

LOG = logging.getLogger(__name__)


###################### SAVE REPEAT #######################################

def save_repeat_fasta(tandem_repeats, file):
    """ save multiple <tandem_repeats> in Fasta format in specified <file>

        At current, only one TR per sequence can be defined, as the identifiers
        in the dict <tandem_repeats> must be unique.

        Parameters: Dict of tandem repeats and identifiers.
            e.g. {'ENSP00012': msa1, 'ENSP00013': msa2}

            >ID
            GHKI
            GHKI
            GH--
    """

    with open(file, 'w', newline='\n') as f:
        for identifier, msa in tandem_repeats.items():
            f.write(">{0}\n".format(identifier))
            f.write("\n".join(msa) + "\n\n")


def save_repeat_stockholm(tandem_repeat, file):
    """ save <tandem_repeat> in STOCKHOLM format in specified <file>

        Parameters: Tandem repeat MSA.
            e.g. ["ACDEF-", "ACCDEF"]

            More information in
            ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf
            STOCKHOLM Format Example:

            # STOCKHOLM 1.0
            seq1 ACDEF...GHIKL
            seq2 ACDEF...GHIKL
            seq3 ...EFMNRGHIKL

            seq1 MNPQTVWY
            seq2 MNPQTVWY
            seq3 MNPQT...
            //

    """

    with open(file, 'w', newline='\n') as fh:
        fh.write("# STOCKHOLM 1.0\n")
        for i, iMSA in enumerate(tandem_repeat):
            fh.write("{} {}\n".format(str(i), iMSA))
        fh.write("//")


def save_repeat_treks(tandem_repeats, file):
    """ Save multiple `tandem_repeats` in T-REKS format in specified `file`

        At current, only one TR per sequence can be defined, as the identifiers
        in the dict <tandem_repeats> must be unique.

        Parameters: Dict of tandem repeats and identifiers.
            e.g. {'ENSP00012': [msa1, begin1], 'ENSP00013': [msa2, begin2]}

            T-REKS format example:

            >a
            Length: 3 residues - nb: 3  from  1 to 10 - Psim:0.8076923076923077 region Length:42
            GHKI
            GHKI
            GH--
            **********************

            Length: 3 residues - nb: 2  from  21 to 27 - Psim:0.7857142857142857 region Length:22
            GHKI
            GH--
            **********************

            >b
            Length: 3 residues - nb: 3  from  1 to 10 - Psim:0.8095238095238095 region Length:20
            UIR
            UIR
            UIR
    """

    with open(file, 'w', newline='\n') as fh:
        for identifier, info in tandem_repeats.items():
            fh.write(">{0}\n".format(identifier))
            fh.write("Length: from {0} to\n".format(str(info[1])))
            fh.write("\n".join(info[0]) + "\n" + "*" * 22 + "\n\n")


############################## READ REPEAT SEQUENCE ######################


def read_fasta(seq_filename, sequence_type='AA'):
    """ Read repeat from file in fasta format.

    Read repeat from file in fasta format.

    Args:
        seq_filename (str):  Path to the repeats containing fasta file
        sequence_type (str): Either "AA" or "DNA"

    Returns:
        (list of Repeat): A list of Repeat instances.
    """

    pat_start = re.compile(r">(.*)")
    pat_repeat_unit = re.compile(r"([\w\.\-]+)")

    # Our possible parser states:
    #
    # 1: searching for sequence name
    # 2: searching for repeat units

    repeats = {}
    state = 1
    with open(seq_filename, "rt") as infile:

        for i, line in enumerate(infile):
            LOG.debug("Line {0}: {1}".format(i, line[0:-1]))
            if 1 == state:
                match = pat_start.match(line)
                if match:
                    LOG.debug(" * (1->2) Found start")
                    LOG.debug("Start: %s", match.group(1))
                    name = match.group(1)
                    repeats[name] = []
                    state = 2

            elif 2 == state:
                match = pat_repeat_unit.match(line)
                if match:
                    LOG.debug(" * (2->2) Found Repeat unit")
                    LOG.debug("Repeat Unit: %s", match.group(1))
                    repeat_unit = match.group(1)
                    repeats[name].append(repeat_unit.replace(".", "-").upper())
                else:
                    LOG.debug(" * (2->1) Found NO further Repeat unit")
                    state = 1

    for i_name, iR in repeats.items():
        yield iR, sequence_type


################################### SIMULATE SEQUENCE ####################

def evolved_tandem_repeats(l, n, n_samples, sequence_type, job_id='job_id',
                           mutation_rate=50, tree='star',
                           indel_rate_per_site=False,
                           return_type='repeat'):
    """ Simulate evolved sequences with ALF.

    Simulate evolved sequences with ALF:
    Dalquen, D. A., Anisimova, M., Gonnet, G. H. & Dessimoz, C.
    ALF--a simulation framework for genome evolution. Molecular Biology and Evolution 29,
    1115–1123 (2012).

    Args:
        l (int): The length of the repeat unit
        n (int): The number of repeat units in the tandem repeat
        n_samples (int):  The number of samples
        sequence_type (str): Either "AA" or "DNA"
        job_id (str): A tag for files produces with ALF, and result files.
        mutation_rate (float): The mutation rate.
        tree (str): The type of tree, e.g. "star" or "birthdeath"
        indel_rate_per_site (int or False): The indel rate per site.
        return_type (str): Either "repeat" or "list"

        sequence_length (int): The total length of the simulated sequence

    Returns:
        Return type depends on ``return_type``: ``Repeat`` or ``Bio.Seq.Seq`` instance.
    """

    runfile_template = os.path.join(DATA_DIR, "ALF", "template.drw")
    alf_exec = os.path.join(EXEC_DIR, "alfsim")
    # create temporary directory
    working_dir = tempfile.mkdtemp()
    LOG.debug("evolvedTR: Created tempfile: %s", working_dir)

    # create working dir
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    # Copy template file and append job specific info
    runfilename = "alf.drw"
    shutil.copyfile(runfile_template, os.path.join(working_dir, runfilename))

    with open(os.path.join(working_dir, runfilename), "a") as runfile:
        if indel_rate_per_site:
            # insertion rate per site and PAM. E.g. for PAM=40 expect
            # aaGainRate*40 insertions per site.
            runfile.write(
                "aaGainRate := " + str(indel_rate_per_site / mutation_rate) + ";\n")
            runfile.write(
                "aaLossRate := " + str(indel_rate_per_site / mutation_rate) + ";\n")
            runfile.write("maxIndelLength := 50;\n")
            runfile.write("indelModel := 'ZIPF';\n")
            runfile.write("Z_c := 1.821;\n")
            runfile.write("DawgPlacement := true;\n")
        runfile.write("uuid := '" + job_id + "';\n")
        runfile.write("mname := '" + job_id + "';\n")
        runfile.write("protStart := " + str(n_samples) + ";\n")
        runfile.write("NSpecies := " + str(n) + ";\n")
        runfile.write("minGeneLength := " + str(l) + ";\n")
        runfile.write("mutRate := " + str(mutation_rate) + ";\n")
        runfile.write("wdir := '" + working_dir + "';\n")

        tree_length = n * mutation_rate
        # parameters concerning the species tree
        if tree == 'star':
            # BDTree, ToLSample, Custom
            runfile.write("treeType := 'Custom':\n")
            runfile.write(
                "tree := MakeStarTree(BirthDeathTree(0.1,0.1," + str(n) + ",10),mutRate):\n")
            runfile.write("treeLength:= %d ;\n" % tree_length)
        else:
            # BDTree, ToLSample, Custom
            runfile.write("treeType := 'BDTree':\n")
            # From the last discussion with Daniel, only the ration of birth to
            # death rate matters, if we scale tree to match pam distance
            runfile.write("birthRate := 0.01:\n")
            runfile.write("deathRate := 0.01:\n")
            # for BDTree: should resulting tree be ultrametric, e.g. all leaves
            # have same distance to origin?
            runfile.write("ultrametric := false:\n")
            runfile.write("treeLength:= %d ;\n" % tree_length)
            # DANIEL: Is a tree := missing?
        if tree not in {'star', 'birthdeath'}:
            LOG.warning(
                "evolvedTR: tree input %s not known, assuming birthdeath tree",
                tree)

        # parameters concerning the substitution models
        if sequence_type == 'AA':
            runfile.write(
                "substModels := [SubstitutionModel('CustomP', ['" +
                os.path.join(
                    DATA_DIR,
                    'ALF',
                    'lg.dat') +
                "'])];\n")
            # CHECK! DOES
            # substModels := [SubstitutionModel('LG')];
            # WORK?
        elif sequence_type == 'DNA':
            runfile.write(
                "substModels := [SubstitutionModel('TN93', [.3, .4, .7],[seq(0.25,4)], true)]:\n")
            # CHECK! DO THE EMPIRICAL PARAMETERS MAKE SENSE AT ALL?
            # compare to http://people.inf.ethz.ch/ddalquen/alf/ALF_manual.pdf
        else:  # CODONS
            runfile.write("substModels := [SubstitutionModel('CPAM')];\n")

    # Determine ALF MSA output file path
    if sequence_type == 'AA':
        infile = os.path.join(working_dir, job_id, "MSA", "MSA_all_aa.phy")
    elif sequence_type == 'DNA':
        infile = os.path.join(working_dir, job_id, "MSA", "MSA_all_dna.phy")
    else:  # CODONS
        infile = os.path.join(
            working_dir,
            job_id,
            "MSA",
            "MSA_all_codon.phy")

    for i in range(10):
        with open(os.path.join(working_dir, "out.txt"), "w") as outfile:
            alf_process = subprocess.Popen([alf_exec,
                                            runfilename],
                                           cwd=working_dir,
                                           stdout=outfile,
                                           stderr=outfile,
                                           close_fds=True)
            alf_process.wait()
        if os.path.isfile(infile):
            break
    else:
        LOG.error('ALFSIM was not able to produce simulated sequence.')

    # shutil.rmtree('/cluster/home/infk/eschaper/spielwiese/')
    #shutil.copytree(working_dir, '/cluster/home/infk/eschaper/spielwiese/')

    """ Read in MSA data from ALF
        The following is a short parser for ALFsim output files
        yielding MSAs as lists of strings,
        based on a special flavour of  Felstein stein MSA files """

    # find a repeat uni
    pattern_start = re.compile("\d+ \d+")
    pattern_seq = re.compile("\S+[ ]+([A-Z\-]+)")

    # Our possible parser states:
    #
    # state 1: Find beginning of MSA (Felsenstein)
    # state 2: Find all repeat units

    # protein ::=
    #    \d \d
    #    \S/\S    repeatunits
    #

    state = 1
    with open(infile, "r") as infile:
        for i, line in enumerate(infile):
            LOG.debug("Line %d: %s", i, line[0:-1])

            if 1 == state:  # Find first repeat unit & save begin
                search = pattern_start.search(line)
                if search:
                    LOG.debug(" *(1->2) Found MSA start")
                    state = 2
                    msa = []

            elif 2 == state:  # Find all repeat units
                search = pattern_seq.search(line)
                if search:
                    LOG.debug(" *(2->2) Found another repeat unit")
                    msa.append(search.group(1))
                else:
                    LOG.debug(" *(2->1) repeat region finished, yielding.")
                    state = 1
                    # YIELD IF WE HAVE FOUND AT LEAST TWO REPEAT UNITS:
                    if len(msa) > 1:
                        if return_type == 'repeat':
                            yield repeat.Repeat(begin=0, msa=msa, sequence_type=sequence_type)
                        elif return_type == 'list':
                            yield msa
                        else:
                            LOG.debug("YIELD: %s",
                                      "".join(msa).replace('-', ''))
                            yield Bio.Seq.Seq("".join(msa).replace('-', ''), sequence_type)

    # delete temporary directory
    try:
        shutil.rmtree(working_dir)
    except OSError:
        pass


def random_sequence(n_samples, sequence_type='AA', return_type='repeat',
                    equilibrium_frequencies='human', l=0, n=0,
                    sequence_length=0):
    """ Simulate random sequence locally.

    Simulate random sequence locally.

    Args:
        n_samples (int):  The number of samples
        sequence_type (str): Either "AA" or "DNA"
        return_type (str): Either "repeat" or "list"
        equilibrium_frequencies (str): Only "human" option available at current
        l (int): The length of the repeat unit
        n (int): The number of repeat units in the tandem repeat
        sequence_length (int): The total length of the simulated sequence

    Returns:
        Return type depends on ``return_type``.
    """

    if sequence_length == 0 and (l == 0 or n == 0):
        LOG.error('The specified sequence_length or the product of l and n was'
                  ' set to 0 for random_sequence simulation')
    else:
        if sequence_length == 0:
            sequence_length = l * n

        if equilibrium_frequencies == 'human':
            file = os.path.join(DATA_DIR, 'Random', "_".join(
                [sequence_type, equilibrium_frequencies, '3']) + '.txt')

        with open(file, 'r') as f:
            a = f.readline()[:-1]
        b = [i.split(' ') for i in a.split('  ')]
        frequencies = {i[0]: int(i[1]) for i in b}
        alphabet = np.unique(i[0] for i in frequencies.keys())

        for _ in range(n_samples):
            seed_int = random.randint(1, sum(frequencies.values()))
            for key, value in frequencies.items():
                seed_int -= value
                if seed_int <= 0:
                    seed = key
                    break

            dimer = [''.join(i) for i in itertools.product(alphabet, repeat=2)]
            third_letter_frequencies = {
                iD: {
                    iA: frequencies[
                        iD +
                        iA] for iA in alphabet} for iD in dimer}

            sequence = seed
            for _ in range(sequence_length - 3):
                next_int = random.randint(
                    1, sum(third_letter_frequencies[sequence[-2:]].values()))
                for key, value in third_letter_frequencies[
                        sequence[-2:]].items():
                    next_int -= value
                    if next_int <= 0:
                        sequence += key
                        break

            # return Repeat instances
            if return_type == 'repeat' and not l == 0 and not n == 0:
                yield [sequence[i * l:(i + 1) * l] for i in range(n)], sequence_type
            elif return_type == 'list':
                yield sequence
            else:  # return seqIO instances
                yield Bio.Seq.Seq(sequence, sequence_type)
