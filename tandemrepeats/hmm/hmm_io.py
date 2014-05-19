import logging
import os, re

logger = logging.getLogger(__name__)

################################## READ HMMMER3  #########################################


def read(hmm_filename, id = None, type = "ACC"):

    '''Read HMM file in HMMER3 format.
       Compare ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf // 8.File formats

       If <id> is defined, return only HMM with an accession tag that matches <id>. <type> identifies the type of accession.
        Store probabilities as log10arithms. [???]

        Parameters:
            <hmm_file> = 'path/to/file.hmm'
            <id> = 'PF00560'


       The general format is::

        HMMER3/b [3.0b2 | June 2009]
        NAME fn3
        ACC PF00041.12
        DESC Fibronectin type III domain
        LENG 86
        ALPH amino
        (...)

        HMM     A       C       D       E       F       G       H       I       (...)       Y
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


       .. important:: The probabilities are stored as negative natural log probabilities. E.g. -ln(0.25) = 1.38629

     '''

    pat_start_HMMER3 = re.compile(r"HMMER3")
    pat_accession = re.compile(r"{}\s+([\w\.]+)".format(type))
    pat_name = re.compile(r"NAME\s+([\w\.]+)")
    pat_HMM = re.compile(r"HMM")
    pat_letters = re.compile(r"\w")
    pat_emissions = re.compile(r"[\w\.]+")
    pat_insertions = re.compile(r"[\d\.]+")
    pat_transition = re.compile(r"[\d\.\*]+")
    pat_end_HMM = re.compile(r"//")

    # Our possible parser states:
    #
    # 0: searching for HMMER3 Tag
    # 0.1 : searching for ACC Tag
    # 1: searching for HMM Tag
    # 2: Skipping line "m->m    m->i    m->d    i->m    i->i    d->m    d->d"
    # 3: searching for emission probabilities OR end of HMM
    # 4: searching for insertion emission probabilities
    # 5: searching for transition probabilities

    lHMM = []

    state = 0
    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            logging.debug("Line {0}: {1}".format(i, line[0:-1]))
            if 0 == state:
                match = pat_start_HMMER3.match(line)
                if match:
                    if id:
                        logging.debug(" * (0->0.1) Found start")
                        state = 0.1
                    else:
                        logging.debug(" * (0->1) Found start")
                        state = 1
            elif 0.1 == state:
                match = pat_accession.match(line)
                if match:
                    iID = match.group(1)
                    if id in iID:
                        logging.debug(" * (0.1->1) Found matching accession.")
                        state = 1
                    else:
                        logging.debug(" * (0.1->0) Found accession. Did not match lID.")
                        state = 0
            elif 1 == state:
                match = pat_HMM.match(line)
                if match:
                    letters = pat_letters.findall(line[3:])
                    logging.debug(" * (1->2) Found HMM start")
                    logging.debug("Letters: {0}".format(letters))
                    hmm = {'letters': letters}
                    state = 2
            elif 2 == state:
                logging.debug(" * (2->3) Skipping line")
                state = 3

            elif 3 == state:
                findall = pat_emissions.findall(line)
                if findall:
                    current_hmm_state = findall[0]
                    if current_hmm_state == 'COMPO':
                        emissions = findall[1:]
                    else:
                        emissions = findall[1:-1]
                    logging.debug(" * (3->4) Found emission probabilities")
                    logging.debug("Current HMM state: {0}".format(current_hmm_state))
                    logging.debug("Emission probabilities: {0}".format(emissions))
                    hmm[current_hmm_state] = {'emissions': emissions}
                    state = 4
                elif pat_end_HMM.match(line):
                    if id:
                        logging.debug(" * (3->TERMINAL) HMM Found and compiled, return HMM.")
                        return hmm
                    else:
                        logging.debug(" * (3->0) Finished HMM compilation")
                        lHMM.append(hmm)
                        state = 0
                else:
                    logging.debug(" * (3->0) Error: No emission line")
                    state = 0

            elif 4 == state:
                findall = pat_insertions.findall(line)
                if findall:
                    logging.debug(" * (4->5) Found insertion emission probabilities")
                    logging.debug("Insertion emission probabilities: {0}".format(findall))
                    hmm[current_hmm_state]['insertion_emissions'] = findall
                    state = 5
                else:
                    logging.debug(" * (4->0) Error: No insertion emission line")
                    state = 0

            elif 5 == state:
                findall = pat_transition.findall(line)
                if findall:
                    logging.debug(" * (5->3) Found transition probabilities")
                    logging.debug("Transition probabilities: {0}".format(findall))
                    hmm[current_hmm_state]['transition'] = findall
                    state = 3
                else:
                    logging.debug(" * (5->0) Error: No transition line")
                    state = 0

    logging.debug(lHMM)
    return lHMM

def split_HMMER_file(hmm_filename, resultdir):

    pat_start_HMMER3 = re.compile(r"HMMER3")
    pat_accession = re.compile(r"ACC\s+([\w\.]+)")
    tmp_file = os.path.join(resultdir, "tmp.hmm")
    state = 0
    fh = None

    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            #logging.debug("Line {0}: {1}".format(i, line[0:-1]))
            if 0 == state:
                match = pat_start_HMMER3.match(line)
                if match:
                    if fh:
                        fh.close()
                        os.rename(tmp_file, os.path.join(resultdir, acc+".hmm"))
                    logging.debug(" * (0->1) Found HMM start")
                    state = 1
                    fh = open(tmp_file, "w")
                fh.write(line)
            elif 1 == state:
                fh.write(line)
                match = pat_accession.match(line)
                if match:
                    acc = match.group(1).split(".")[0]
                    logging.debug(" * (1->0) Found accession")
                    state = 0

def read_HMMER_acc_lengths(hmm_filename):

    '''Read HMM file in HMMER3 format. (See definition further down)
        Return the PFAM ID and the lengths of each model.
    '''

    pat_accession = re.compile(r"ACC\s+([\w\.]+)")
    pat_length = re.compile(r"LENG\s+([\w\.]+)")
    state = 0
    d = {}

    with open(hmm_filename, "rt") as infile:

        for i, line in enumerate(infile):
            #logging.debug("Line {0}: {1}".format(i, line[0:-1]))
            if 0 == state:
                match = pat_accession.match(line)
                if match:
                    acc = match.group(1)
                    logging.debug(" * (0->1) Found name")
                    state = 1
            elif 1 == state:
                match = pat_length.match(line)
                if match:
                    d[acc] = int(match.group(1))
                    logging.debug(" * (1->0) Found length")
                    state = 0

    return d

# def write_HMMER(hmm, hmm_filename):
#
#     ''' Write <hmm> too <hmm_filename> in HMMER3 format.
#        Compare ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf // 8.File formats '''
#
#
#
#
# def read_HMM_smart(hmm_filename, id = None, type = "ACC"):
#
#     ''' Write <hmm> too <hmm_filename> in smart format. '''
#     # Format not decided yet.


##################################### Main ###############################################

def main():
    lHMM = read_HMMER('you/favourite/path')

if __name__ == "__main__":
    main()
