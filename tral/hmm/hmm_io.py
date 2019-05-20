# (C) 2015 Elke Schaper

"""
:synopsis: Input/output for hmms.

.. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import os
import re
import math
import tral

LOG = logging.getLogger(__name__)

# ################################ READ HMMMER3  #########################


def read(hmm_filename, id=None):
    """Read HMM file in HMMER3 format.

    HMMER3 file format is described in detail in
    ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf,
    section 8.

    The general file format is::

        HMMER3/b [3.0b2 | June 2009]
        NAME fn3
        ACC PF00041.12
        DESC Fibronectin type III domain
        LENG 86
        ALPH amino
        (...)

        HMM          A       C       D       E       F       G       H       I    (...)    Y
                    m->m    m->i    m->d    i->m    i->i    d->m    d->d
          COMPO   2.70271 4.89246 3.02314 2.64362 3.59817 2.82566 3.74147 3.08574 (...) 3.22607
                  2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                  0.00338 6.08833 6.81068 0.61958 0.77255 0.00000       *
          1       3.16986 5.21447 4.52134 3.29953 4.34285 4.18764 4.30886 3.35801 (...) 3.93889 1 - -
                  2.68629 4.42236 2.77530 2.73088 3.46365 2.40512 3.72505 3.29365 (...) 3.61514
                  0.09796 2.38361 6.81068 0.10064 2.34607 0.48576 0.95510
          2       2.70230 5.97353 2.24744 2.62947 5.31433 2.60356 4.43584 4.79731 (...) 4.25623 3 - -
                  2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                  0.00338 6.08833 6.81068 0.61958 0.77255 0.48576 0.95510
          (...)
          86      3.03720 5.94099 3.75455 2.96917 5.26587 2.91682 3.66571 4.11840 (...) 4.99111 121 - E
                  2.68618 4.42225 2.77519 2.73123 3.46354 2.40513 3.72494 3.29354 (...) 3.61503
                  0.00227 6.08723       * 0.61958 0.77255 0.00000       *
        //

    .. important:: The probabilities in this format are stored as negative
        natural log probabilities, e.g. -ln(0.25) = 1.38629. The special case
        of 0 probability is stored as ``*`` which in fact means -∞ (minus
        infinity).

    Args:
        hmm_filename (str): Path to the file with model data in the HMMER3
            file format.
        id (str, optional): The identifier for the model to be returned.
            E.g. for Pfam the ``id`` can look like this: 'PF00560'.
            If defined, the function returns only the HMM with an identifier
            that matches the provided id. If a model in the file does not have
            an identifier to check against, it is skipped.

    Returns:
        dict: A dictionary of parameters required to initialise the HMM
        including the id, alphabet, emission probabilities and transition
        probabilities.

        Output format::

            {
                'id': 'PF08261.7',
                'letters': ['A', 'C', 'D', ..., 'Y'],
                'COMPO':
                    {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                      'emissions': [2.28205, 5.14899, ..., 1.92022],
                      'transitions': [0.01467, 4.62483, ..., -inf]},
                '1':
                    {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                     'emissions': [1.00089, 4.54999, ..., 5.23581],
                     'transitions': [0.01467, 4.62483, ..., 0.95510]},
                ...
                '8':
                    {'insertion_emissions': [2.68618, 4.42225, ..., 3.61503],
                     'emissions': [4.12723, 5.39816, ..., 4.58094],
                     'transitions': [0.00990, 4.62006, ..., -inf]}
            }

    """
    pat_start_HMMER3 = re.compile(r"HMMER3")
    pat_accession = re.compile(r"ACC\s+([\w\.]+)")
    pat_HMM = re.compile(r"HMM")
    pat_letters = re.compile(r"\w")
    pat_emissions = re.compile(r"[\w\.]+")
    pat_insertions = re.compile(r"[\d\.]+")
    pat_transition = re.compile(r"[\d\.\*]+")
    pat_end = re.compile(r"//")

    # Our possible parser states:
    #
    # 0: searching for HMMER3 Tag
    # 0.1 : searching for identifier Tag
    # 1: searching for HMM Tag
    # 2: Skipping line "m->m    m->i    m->d    i->m    i->i    d->m    d->d"
    # 3: searching for emission probabilities OR end of HMM
    # 4: searching for insertion emission probabilities
    # 5: searching for transition probabilities

    state = 0
    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            LOG.debug("Line %s: %s", i, line[0:-1])
            if 0 == state:
                match = pat_start_HMMER3.match(line)
                if match:
                    LOG.debug(" * (0->0.1) Found model start.")
                    state = 0.1

            elif 0.1 == state:
                hmm = {'id': None}
                match = pat_accession.match(line)
                if match:
                    iID = match.group(1)
                    if id:
                        if id in iID:
                            LOG.debug(" * (0.1->1) Found matching"
                                      " identifier.")
                            hmm['id'] = iID
                            state = 1
                        else:
                            LOG.debug(" * (0.1->0) Found identifier which "
                                      "does not the one provided, skipping "
                                      "model.")
                            state = 0
                    else:
                        hmm['id'] = iID
                        state = 1

                # This is done in order to be able to parse HMMs without an
                # accession number, e.g. when created from a Repeat.
                match2 = pat_HMM.match(line)
                if match2:
                    letters = pat_letters.findall(line[3:])
                    if id:
                        LOG.debug(" * (0.1->0) No identifier found to "
                                  "compare to, skipping model.")
                        state = 0
                    else:
                        LOG.debug(" * (0.1->2) No identifier, found"
                                  " HMM start")
                        LOG.debug("Letters: %s", letters)
                        hmm['letters'] = letters
                        state = 2

            elif 1 == state:
                match = pat_HMM.match(line)
                if match:
                    letters = pat_letters.findall(line[3:])
                    LOG.debug(" * (1->2) Found HMM start")
                    LOG.debug("Letters: %s", letters)
                    hmm['letters'] = letters
                    state = 2

            elif 2 == state:
                LOG.debug(" * (2->3) Skipping line")
                state = 3

            elif 3 == state:
                findall = pat_emissions.findall(line)
                if findall:
                    current_hmm_state = findall[0]
                    if current_hmm_state == 'COMPO':
                        string_emissions = findall[1:]
                    else:
                        string_emissions = findall[1:1 + len(hmm['letters'])]
                    LOG.debug(" * (3->4) Found emission probabilities")
                    LOG.debug("Current HMM state: %s", current_hmm_state)
                    emissions = [float(i) if i != '*' else -float('inf')
                                 for i in string_emissions]

                    LOG.debug("Emission probabilities: %s", emissions)
                    hmm[current_hmm_state] = {'emissions': emissions}
                    state = 4
                elif pat_end.match(line):
                    if id:
                        LOG.debug(" * (3->TERMINAL) HMM Found and compiled,"
                                  " return HMM.")
                        LOG.info("Yielding {}, stopping".format(hmm['id']))
                        yield hmm
                    else:
                        LOG.debug(" * (3->0) Finished HMM compilation")
                        LOG.info("Yielding {}".format(hmm['id']))
                        yield hmm
                        state = 0
                else:
                    LOG.debug(" * (3->0) Error: No emission line")
                    state = 0

            elif 4 == state:
                findall = pat_insertions.findall(line)
                if findall:
                    LOG.debug(" * (4->5) Found insertion emission"
                              " probabilities")
                    emissions = [float(i) if i != '*' else -float('inf')
                                 for i in findall]
                    LOG.debug("Insertion emission probabilities: %s",
                              emissions)
                    hmm[current_hmm_state]['insertion_emissions'] = emissions
                    state = 5
                else:
                    LOG.debug(" * (4->0) Error: No insertion emission line")
                    state = 0

            elif 5 == state:
                findall = pat_transition.findall(line)
                if findall:
                    LOG.debug(" * (5->3) Found transition probabilities")
                    transitions = [float(i) if i != '*' else -float('inf')
                                   for i in findall]
                    LOG.debug("Transition probabilities: %s", transitions)
                    hmm[current_hmm_state]['transition'] = transitions
                    state = 3
                else:
                    LOG.debug(" * (5->0) Error: No transition line")
                    state = 0


def split_HMMER3_file(hmm_filename, resultdir):
    """Split HMMER3 models from a single file ``hmm_filename`` into many files
    in ``resultdir``.

    Helper function: split HMMER3 models from a single file ``hmm_filename``
    into many files in ``resultdir``. The models are named after the HMM
    accession.

    Args:
      hmm_filename (str): Path to HMMER3 file.
      resultdir (str): Path to directory where result files are stored.

    """
    pat_start_HMMER3 = re.compile(r"HMMER3")
    pat_accession = re.compile(r"ACC\s+([\w\.]+)")
    tmp_file = os.path.join(resultdir, "tmp.hmm")
    state = 0
    fh = None
    acc = None

    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            LOG.debug("Line %s: %s", i, line[0:-1])
            if 0 == state:
                match = pat_start_HMMER3.match(line)
                if match:
                    if fh:
                        fh.close()
                        os.rename(
                            tmp_file,
                            os.path.join(
                                resultdir,
                                acc +
                                ".hmm"))
                    LOG.debug(" * (0->1) Found HMM start")
                    state = 1
                    fh = open(tmp_file, "w")
                fh.write(line)
            elif 1 == state:
                fh.write(line)
                match = pat_accession.match(line)
                if match:
                    acc = match.group(1).split(".")[0]
                    LOG.debug(" * (1->0) Found accession")
                    state = 0


def read_HMMER_acc_lengths(hmm_filename):
    """Read HMM file in HMMER3 format. Return the PFAM ID and the lengths of
    each model.

    Read HMM file in HMMER3 format. Return the PFAM ID and the lengths of
    each model.

    Args:
      hmm_filename (str): Path to HMMER3 file.

    Returns:
        (dict of str: int): The number of match states for all HMM models in
        ``hmm_filename``.

      ..  todo:: Decide whether this function is needed.

    """
    pat_accession = re.compile(r"ACC\s+([\w\.]+)")
    pat_length = re.compile(r"LENG\s+([\w\.]+)")
    state = 0
    data = {}

    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            LOG.debug("Line %s: %s", i, line[0:-1])
            if 0 == state:
                match = pat_accession.match(line)
                if match:
                    acc = match.group(1)
                    LOG.debug(" * (0->1) Found name")
                    state = 1
            elif 1 == state:
                match = pat_length.match(line)
                if match:
                    data[acc] = int(match.group(1))
                    LOG.debug(" * (1->0) Found length")
                    state = 0

    return data


def to_fixed_width(n, max_width, allow_overflow=True, do_round=True):
    """Format a float to a fixed-width string.

    Args:
      - n (float):              number to format
      - max_width (int):        width of output string
      - allow_overflow (bool):  allow the result to be larger than max_width for
                                large n? Raises OverflowError if False
      - do_round (bool):        Round final digit?
    Returns:
        (str)   string representation of n, either padded or rounded to have
                max_width characters
    """
    # Adapted from https://stackoverflow.com/a/43397325/81658
    if do_round:
        for i in range(max_width - 2, -1, -1):
            str0 = '{:.{}f}'.format(n, i)
            if len(str0) <= max_width:
                break
    else:
        str0 = '{:.42f}'.format(n)
        int_part_len = str0.index('.')
        if int_part_len <= max_width - 2:
            str0 = str0[:max_width]
        else:
            str0 = str0[:int_part_len]
    if (not allow_overflow) and (len(str0) > max_width):
        raise OverflowError("Impossible to represent in fixed-width non-scientific format")
    return str0


def write_HMMER(hmm, hmm_file):
    """Write <hmm> too <hmm_filename> in HMMER3 format.

    See Section 8 (File Formats) of
    ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf
    for the file format specification.

    Args:
        hmm (HMM): HMM object to write
        hmm_file (str or file): filename or file-like object to write to
    """
    close = False
    out = None
    try:
        # file or filename
        if hasattr(hmm_file, 'write'):
            out = hmm_file
        else:
            close = True
            out = open(hmm_file, 'w')

        # Header section
        out.write("HMMER3/f [TRAL {}]\n".format(tral.__version__))

        if hmm.id:
            name = hmm.id
        elif type(hmm_file) == str:
            # use filename
            name = os.path.splitext(os.path.basename(hmm_file))[0]
        elif hasattr(hmm_file, "name"):
            name = os.path.splitext(os.path.basename(hmm_file.name))[0]
        else:
            # final fallback
            name = "TRAL"
        out.write("NAME  {}\n".format(name))

        if hmm.id:
            out.write("ACC   {}\n".format(hmm.id))

        out.write("LENG  {}\n".format(hmm.l_effective))

        if len(hmm.alphabet) == 20:
            alphabet = "amino"
        elif len(hmm.alphabet) == 4:
            if set(hmm.alphabet) == {"A", "G", "C", "T"}:
                alphabet = "DNA"
            elif set(hmm.alphabet) == {"A", "G", "C", "U"}:
                alphabet = "RNA"
            else:
                alphabet = "custom"
        else:
            alphabet = "custom"

        out.write("ALPH  {}\n".format(alphabet))
        out.write("MAP   yes\n")
        out.write("CONS  no\n")
        out.write("RF    no\n")
        out.write("MM    no\n")
        out.write("CS    no\n")
        # HMM header line
        out.write("HMM       ")
        out.write(" ".join([l.center(7) for l in hmm.alphabet]))
        out.write("\n")
        out.write(" " * 10)
        transitions = ["m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d"]
        out.write(" ".join([l.center(7) for l in transitions]))
        out.write("\n")

        # Model section
        # TODO is hmm.hmmer always consistent with hmm.p_0/p_e/p_t?
        if not hmm.hmmer:
            raise Exception("Unsupported: HMM lacks hmmer-style states")

        # Composition
        if "COMPO" not in hmm.hmmer:
            # TODO Read BEGIN probabilites independently of COMPO line
            raise KeyError("COMPO currently required")

        def fmt_prob(p):
            if p > 1e8 or math.isinf(p):
                return "      *"
            else:
                return to_fixed_width(p, 7)

        out.write("  COMPO   ")
        out.write("  ".join([fmt_prob(p) for p in hmm.hmmer['COMPO']['emissions']]))
        out.write("\n")
        # BEGIN probabilities
        out.write(" " * 10)
        out.write("  ".join([fmt_prob(p) for p in hmm.hmmer['COMPO']['insertion_emissions']]))
        out.write("\n")
        out.write(" " * 10)
        out.write("  ".join([fmt_prob(p) for p in hmm.hmmer['COMPO']['transition']]))
        out.write("\n")

        # Main states
        for i in range(1, hmm.l_effective + 1):
            out.write("{:>7}   ".format(i))
            out.write("  ".join([fmt_prob(p) for p in hmm.hmmer[str(i)]['emissions']]))
            out.write(" {map} {cons} {rf} {mm} {cs}\n".format(map=i, cons="-", rf="-", mm="-", cs="-"))
            out.write(" " * 10)
            out.write("  ".join([fmt_prob(p) for p in hmm.hmmer[str(i)]['insertion_emissions']]))
            out.write("\n")
            out.write(" " * 10)
            out.write("  ".join([fmt_prob(p) for p in hmm.hmmer[str(i)]['transition']]))
            out.write("\n")
        out.write("//")

    finally:
        if close and out:
            out.close()


#
#
# def read_HMM_smart(hmm_filename, id = None, type = "ACC"):
#
#     ''' Write <hmm> too <hmm_filename> in smart format. '''
#     # Format not decided yet.
