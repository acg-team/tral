# (C) 2011, Alexander Korsunsky
# (C) 2011-2015 Elke Schaper

"""
    :synopsis: Execution of repeat detection algorithms

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

from collections import OrderedDict
import distutils
import itertools
import logging
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile

from tral import configuration
from tral.sequence import repeat_detection_io

LOG = logging.getLogger(__name__)

CONFIG_GENERAL = configuration.Configuration.instance().config
REPEAT_DETECTOR_PATH = CONFIG_GENERAL["sequence"]["repeat_detector_path"]


class BinaryExecutable:

    """ Contains the executable, and combines executable with parameters.

    Contains the executable, and combines executable with parameters.

     Attributes:
        binary(str): Path to binary

    """

    def __init__(self, binary=None, name=None):
        """ Construct a BinaryExecutable object.

        Construct a BinaryExecutable object.

        Args:
            binary (str): Path to the binary
            name(str): Name of the programme.
        """

        if not binary:
            raise TypeError("A binary executable must be provided.")
        try:
            self.binary = shutil.which(binary)
        except:
            self.binary = distutils.spawn.find_executable(binary)
        if not self.binary:
            raise ValueError(
                "The executable {} does not exist, although {} was selected "
                "to be executed. Please make sure the executable is in the system path, or "
                "the path to the executable is correctly set in config.ini".format(
                    binary, name))

    def get_execute_tokens(self, *args):
        """Return the tokens to invoke the program with the arguments args"""

        return [self.binary] + list(args)

    def get_execute_line(self, *args):
        """Return the command line to invoke the program with the arguments args"""
        return " ".join(self.get_execute_tokens(*args))


def check_java_errors(outfile, errfile, log=LOG, procname=None):
    """ Check for java problems. Return True if there were problems, else False.

    Check for these java errors:

     * Stdout file is empty but stderr file is not.
     * Java Exception string is indicated in the errfile

     Return True if there were problems, else False.

    Args:
        outfile (file handle):  Redirected standard output channel file.
        errfile (file handle): Redirected standard error channel file
            If None, no copies are saved.
        log (?): Name of the log to issue log messages to. If none, no log
                messages will be issued.

    .. todo:: Complete docstring
    """

    has_error = False

    with open(outfile, mode='r') as of, open(errfile, mode='r') as ef:
        outfile_size = of.seek(0, os.SEEK_END)
        errfile_str = ef.read()

    errfile_size = len(errfile_str)

    if errfile_size != 0 and outfile_size == 0:
        has_error = True
        LOG.warning(
            "Process \"%s\" has empty STDOUT but non-empty STDERR",
            procname)

    pat_javaexc = re.compile(
        r'Exception in thread ".+" \S+:.*$(\s+at .+)+',
        re.M)

    m = pat_javaexc.search(errfile_str)
    if m:
        has_error = True
        LOG.warning(
            "Java Exception probably occured in process \"%s\"! Exception information:\n%s",
            procname,
            m.group(0))

    return has_error


class TRDetector(object):

    def __init__(self, executable):
        """Construct a TRDetector object with executable als executable object"""
        LOG.debug(executable)
        self.__executable = executable

    def run_process(self, working_dir, *args):
        """Launch a detector process

        Arguments:
        working_dir -- Working directory to run process in
        args -- Arguments to pass to the process

        Return Value:
        A tuple with:
            - The return code
            - The name to the file with the redirected standard output channel (stdout)
            - The name to the file with the redirected standard error channel (stderror)

        The standard output (stdout) and standard error (stderr) channels will
        be redirected to working_dir/stdout.txt and working_dir/stdout.txt
        respectively.
        """

        # create working directory in top level working directory,
        # if it does not yet exist
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)

        stdoutfname = os.path.join(working_dir, "stdout.txt")
        stderrfname = os.path.join(working_dir, "stderr.txt")

        # open stdout and stderr output file
        __stdout_file = open(stdoutfname, mode='w')
        __stderr_file = open(stderrfname, mode='w')

        LOG.debug("Launching process: %s in %s",
                  self.__executable.get_execute_line(*args), working_dir)

        LOG.debug("Launching process tokens: %s in %s",
                  self.__executable.get_execute_tokens(*args), working_dir)
        # launch process
        __process = subprocess.Popen(
            self.__executable.get_execute_tokens(*args), cwd=working_dir,
            stdout=__stdout_file, stderr=__stderr_file, close_fds=True,
        )
        __process.wait()

        __stdout_file.close()
        __stderr_file.close()

        return __process.returncode, stdoutfname, stderrfname


class DetectorHHrepID(TRDetector):
    name = 'HHrepID'
    displayname = "HHrepID"

    """ Execute vi
    ./hhrepid_32 -i infile -v 0 -nofilt -d dummyHMM.hmm -o resultfile
    ./hhrepid -i <query> -d <path to cal.hhm> -tp <path to tp.dat> -fp <path to fp.dat> """

    class Configuration:

        def __init__(self):
            self.boolopts = {
                # turn off low-complexity filter (default: filter on)
                "-nofilt": False,
                # calibration by shuffling instead of database search (default:
                # off)
                "-shuffle": False
            }

            self.valopts = {
                # <file> input query alignment  (fasta/a2m/a3m) or HMM file (.hhm)
                "-i": None,
                "-d": REPEAT_DETECTOR_PATH['HHrepID_dummyhmm'],   # <path> dummy hmm database file
                "-o": 'hhrepID.o',    # <file> write results and multiple sequence alignment to file (default=none)
                "-v": 0,           # -v: verbose mode (default: show only warnings)  ;  -v 0: suppress all screen outpu
                "-P": None,        # <float> max p-value of suboptimal alignments in all search rounds but the last one (def=0.1)
                "-R": None,        # <float> max p-value of repeats  (def=1)
                "-T": None,        # <float> max total repeat p-value (def=0.001)
                "-alpha": None,        # <float> For calculating repeat p-values: weight of n'th suboptimal alignment vs. (n+1)-st suboptimal alignmen
                "-k": None,        # [0,inf[ For calculating repeat p-values: maximal number of suboptimal alignments considered (def=1)
                "-cont": None,        # <float>  probability threshold for masking of inconsistent cells (def=0.001)
                "-mrgr": None,        # [0,inf[  number of merge rounds for achieving consistency (def=3)
                "-lcon": None,        # <float>  preserved local context in shuffling (def=2)
                "-lmin": None,        # [0,inf[  minimal length of repeats to be identified (def=7)

            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            toks = []
            if infile:
                toks += ['-i', infile]

            for optstring, optvalue in self.valopts.items():
                if optvalue is not None:
                    toks += [str(optstring), str(optvalue)]

            for optstring, optvalue in self.boolopts.items():
                if optvalue:
                    toks.append(optstring)

            return toks

    def __init__(self, name=name):
        """Construct DetectorHHrepID object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorHHrepID, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, DetectorHHrepID.name)

        prog_args = self.config.tokens(infile=infile)

        # execute detector process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args)

        # CHECK WHETHER outfile exists, else give warning and return empty
        # list.

        # Process output file, return results
        resultfile = os.path.join(wd, self.config.valopts["-o"])
        if os.path.isfile(resultfile):
            with open(resultfile, "r") as infilehandle:
                tmp = list(
                    repeat_detection_io.hhpredid_get_repeats(infilehandle))
            return tmp
        else:
            return []


class DetectorPhobos(TRDetector):
    name = 'PHOBOS'
    displayname = "PHOBOS"

    # execute via phobos --printRepeatSeqMode 3 DNA.faa phobos.o
    # or phobos --help to see all options

    class Configuration:

        def __init__(self):
            self.boolopts = {
                # Phobos will choose the shorter of the two alignments instead
                # of the longer. Phobos will run faster with this option.
                "--preferShorterRepeats": False,
                # dontRemoveMostlyOverlapping. Default: DO remove mostly
                # overlapping repeats in favour of the one with the higher
                # internal score
                "--dontRemoveMostlyOverlapping": False,
            }

            self.scriptopts = {
                "outputfile": 'phobos.o',
            }

            self.valopts = {
                # <int> Display of TR. 0: asIs, 1: Alphabetical normal form, 2: Alphabetical normal form also considering the reverse complement. Default: 2
                "--reportUnit": None,
                "--NPerfectionMode": None,  # <int> Treatment of N.  0: asMismatch, 1: asNeutral, 2: asMatch. Default: 0.
                "--printRepeatSeqMode": '3',  # <int> Show the repeat unit alignment
                "--outputFormat": None,  # <int> 0 : Use the Phobos output format. Default: 0
                "--minPerfection": None,  # <float> Minimum perfection of a satellites. Default: 0.
                "--maxPerfection": None,  # <float> Maximum perfection of a satellites. Default: 100.
                "--maximum_score_reduction": None,  # The maximum amount the score can be reduced before search is aborded. Typical: 6*mismatch-penalty or infinite. Default: infinite
                "--recursion": None,  # The recursion depth used in the search. Values in the range 3 to 7 are recommended. A value of 0 implies a search for perfect repeats only. Default: Not stated
                "--minUnitLen": None,  # <int> Minimum unit length. Default: 1
                "--maxUnitLen": None,  # <int> Maximum unit length. Default: 10
                "--minLength_b": None,  # <float> The minimum length of a repeat is determined with: maximum( minLength, minLength_a + minLength_b*(unit-length) ). Default value of minLength_b: 0
                "--minLength_a": None,  # <int> Default value of minLength_a: 0
                "--minLength": None,  # <int> Default value of minLength: 0
                "--minScore_b": None,  # <float> The minimum score of a repeat is determined with: maximum( minScore, minScore_a + minScore_b*(unit-length) ). Default value of minScore_b: 1
                "--minScore_a": None,  # <int> Default value of minScore_a: 0
                "--minScore": None,  # <int> Default value of minScore: 6
                "--mismatchScore": None,  # <int> Score for mismatch - must be negative. Default: -5. Match score is fixed to one.
                "--indelScore": None,  # <int> Score for indels - must be negative. Default: -5. Match score is fixed to one.
                "--searchMode": 'imperfect'  # ['exact','extendExact','imperfect']
            }

        def set_working_dir(self, working_dir):
            self.scriptopts['outputfile'] = os.path.join(
                working_dir,
                self.scriptopts['outputfile'])
            if not os.path.isdir(working_dir):
                os.makedirs(working_dir)

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            toks = [
                optstring for optstring,
                optvalue in self.boolopts.items() if optvalue]

            for optstring, optvalue in self.valopts.items():
                if optvalue is not None:
                    toks.append(optstring)
                    toks.append(str(optvalue))

            if infile:
                toks.append(infile)
            toks.append(self.scriptopts['outputfile'])
            return toks

    def __init__(self, name=name):
        """Construct DetectorPhobos object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorPhobos, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, DetectorPhobos.name)
        self.config.set_working_dir(working_dir=wd)

        # execute detector process
        retcode, stdoutfname, stderrfname = super(
            DetectorPhobos, self).run_process(
            wd, *self.config.tokens(infile))
        # alternatively:
        #prog_args = self.config.tokens(infile=infile)
        #retcode, stdoutfname, stderrfname = super().run_process(wd, *prog_args)

        # Process output file, return results If outfile does not exist return
        # empty list
        if os.path.isfile(self.config.scriptopts['outputfile']):
            with open(self.config.scriptopts['outputfile'], "r") as resultfilehandle:
                tmp = list(
                    repeat_detection_io.phobos_get_repeats(resultfilehandle))
            return tmp
        else:
            LOG.warning(
                "Did not find Phobos result file in %s",
                self.config.scriptopts['outputfile'])
            return []


class DetectorTRED(TRDetector):
    name = 'TRED'
    displayname = "TRED"

    """ tred is a shell script executing:
            tred1 myDNA.faa intermediate_output
            tred2 myDNA.faa intermediate_output output1 output2
        At the moment, we do not know how to add parameters to TRED    """

    class Configuration:

        def __init__(self):
            # For now, use Tred with default options only
            self.boolopts = {
            }

            self.valopts = {
            }

        def set_result_file(self, result_file=None):
            self.result_file = result_file

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            toks = []
            if infile:
                toks = [infile]
            toks.append(self.result_file)

            return toks

    def __init__(self, name=name):
        """Construct DetectorTRED object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorTRED, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, DetectorTRED.name)

        self.result_file = os.path.join(wd, 'tred.o')

        self.config.set_result_file(result_file=self.result_file)
        prog_args = self.config.tokens(infile=infile)

        # execute detector process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args)

        # CHECK WHETHER outfile exists, else give warning and return empty
        # list.

        # Process output file, return results
        if os.path.isfile(self.result_file):
            with open(self.result_file, "r") as infilehandle:
                tmp = list(repeat_detection_io.tred_get_repeats(infilehandle))
            return tmp
        else:
            LOG.warning(
                "Did not find TRED result file in %s",
                self.result_file)
            return []


class DetectorTREKS(TRDetector):
    name = 'T-REKS'
    displayname = "T-REKS"

    class Configuration:

        def __init__(self):
            self.boolopts = {
                "-overlapfilter": True,
                "-nosplit": True
            }

            self.valopts = {
                "-msaMode": None,
                "-clustal": None,
                "-db": None,
                "-muscle": None,
                "-outfile": None,
                "-similarity": None,
                "-kmeans": None,
                "-varIndels": None,
                "-type": None
            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            bool_toks = [
                optstring for optstring,
                optvalue in self.boolopts.items() if optvalue]
            value_toks = ([
                optstring + "=" + str(optvalue)
                for optstring, optvalue in self.valopts.items()
                if optvalue is not None
            ])

            if infile:
                value_toks.append("-infile=" + infile)

            return bool_toks + value_toks

    def __init__(self,
                 name=name,
                 ):
        """Construct DetectorTREKS object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorTREKS, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, DetectorTREKS.name)

        prog_args = self.config.tokens(infile=infile)

        # execute detector process
        retcode, stdoutfname, stderrfname = \
            super(DetectorTREKS, self).run_process(wd, *prog_args)

        if check_java_errors(stdoutfname, stderrfname,
                             log=LOG, procname=DetectorTREKS.displayname):
            return []

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            tmp = list(repeat_detection_io.treks_get_repeats(outfile))
        return tmp


class DetectorTRF(TRDetector):
    name = 'TRF'
    displayname = "TRF_Benson"

    # Execute via:
    # ./trf404.linux64.exe inputfile 2 7 7 80 10 50 500 -d
    # Result file in this case is called inputfile.2.7.7.80.10.50.500.1.txt.html
    # WARNING: This result file name is only used if there is only one sequence in the input file!
    ## result_file_name = ".".join([inputfile] + [str(optvalue) for optstring, optvalue in self.valopts.items() if optvalue != None] + ['1.txt.html'])

    # ./trf404.linux64.exe File Match Mismatch Delta PM PI Minscore MaxPeriod [options]

    class Configuration:

        def __init__(self):
            self.boolopts = {
                # no redundancy elimination # Already without the -r flag set,
                # the three best redundant solutions are displayed
                "-r": False,
                #"-h"            : False,    # When -h is set, the output html of interest is not written :(
                #"-d"            : False,    # produce data file, not needed at curren
                #"-m"            : False,    # masked sequence file
                #"-f"            : False,    # flanking sequence
            }

            # The order of TRF flags matters. Hence, we use an ordered
            # dictionary for valopts.
            valopts_tmp = [
                ("File", None),  # input file
                ("Match", 2),     # matching weigh
                ("Mismatch", 7),     # mismatching penalty
                ("Delta", 7),     # indel penalty
                ("PM", 80),    # integer match probability
                ("PI", 10),    # integer indel probability
                ("Minscore", 50),    # minimum alignment score to report
                ("MaxPeriod", 1000)    # maximum period size to report

            ]
            self.valopts = OrderedDict(valopts_tmp)

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""

            toks = []
            if infile:
                toks.append(infile)

            for optstring, optvalue in self.valopts.items():
                if optvalue is not None:
                    toks.append(str(optvalue))

            for optstring, optvalue in self.boolopts.items():
                if optvalue:
                    toks.append(optstring)

            return toks

    def __init__(self, name=name):
        """Construct DetectorTRF object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorTRF, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, DetectorTRF.name)

        prog_args = self.config.tokens(infile=infile)

        # execute detector process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args)

        # get the infile name from the infile path
        infile_name = re.findall(r"/([_\.\w]+)$", infile)[0]

        # WARNING! TRF is a volatile programme, this result_file_name might not always be true.
        # Maybe a general search on txt.html elements in the result dir might
        # be more appropiate.
        result_file_name = ".".join(
            [infile_name] +
            [
                str(optvalue) for optstring,
                optvalue in self.config.valopts.items() if optvalue is not None] +
            ['1.txt.html'])

        # Process output file, return results
        with open(os.path.join(wd, result_file_name), "r") as outfile:
            tmp = list(repeat_detection_io.trf_get_repeats(outfile))
        return tmp


class DetectorTrust(TRDetector):
    name = 'TRUST'
    displayname = "TRUST"

    class Configuration:

        def __init__(self):
            self.boolopts = {
                "-noseg": True,
                "-force": False,
            }

            self.valopts = {
                "-matrix": REPEAT_DETECTOR_PATH['TRUST_substitutionmatrix'],
                "-gapo": "8",
                "-gapx": "2",
                "-procTotal": "1",
            }  # procTotal: When running on a cluster, specify the amount of processors used  # procNr cpu_number:  When running on a cluster, specify 0-based CPU-number and amount of processors used

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            toks = [
                optstring for optstring,
                optvalue in self.boolopts.items() if optvalue]

            for optstring, optvalue in self.valopts.items():
                if optvalue is not None:
                    toks.append(optstring)
                    toks.append(str(optvalue))

            if infile:
                toks.append("-fasta")
                toks.append(infile)

            return toks

    def __init__(self, name=name):
        """Construct DetectorTrust object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorTrust, self).__init__(executable)

    def run_process(self, working_dir, infile):
        """Run detector process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, DetectorTrust.name)

        # execute detector process
        retcode, stdoutfname, stderrfname = super(
            DetectorTrust, self).run_process(
            wd, *self.config.tokens(infile))

        if check_java_errors(stdoutfname, stderrfname,
                             log=LOG, procname=DetectorTrust.displayname):
            return []

        # shutil.copyfile(stdoutfname, YOUR FAVOURITE PATH)

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            return list(repeat_detection_io.trust_get_repeats(outfile))


class DetectorXStream(TRDetector):
    name = 'XSTREAM'
    displayname = "XSTREAM"

    class Configuration:

        def __init__(self):
            self.boolopts = {
                #"-n" : False,   # sequence type: nucleotide
                "-l": False,   # don't look for periods >1000
                #"-p" : False,   # don't print repeat alignmen
                "-h": True,    # write to console, not HTML
                #"-C" : False,   # don't print any repeat info
                "-z": True,  # generate spreadsheet .xls output
                #"-Z" : False,   # create statistics spreadshee
                #"-s" : False,   # parse output by input sequence
                "-f": False,   # perform divide & conquer on inpu
                #"-G" : False,   # create color-coded TR block diagram
                #"-G*" : False,  # create single color TR block diagram
                #"-B" : False,   # create sequence comparison map PNG
                #"-O" : False,   # print "mismatch/gaps colored" HTML outpu
                "-N": False,   # turn off nesting
                "-o": False,   # don't remove overlapping TR domain
                #"-S" : False,   # input strain into database
                "-V": True,   # don't invoke D & C automatically for long inpu
            }

            self.valopts = {
                "-g": None,    # [integer] set gap number
                "-i": None,    # [0-1] set matching threshold for detection
                "-I": None,    # [0-1] set cons error threshold for printing
                "-D": None,    # [0-1] set max % indels for printing
                "-m": 1,    # [integer] set minimum period
                "-x": None,    # [integer] set maximum period
                "-e": None,    # [any number>=2] set minimum copy number
                # "-a" : None,    # [string] add identifier string to output file
                "-f": None,    # [integer] set fragment length
                #"-d" : None,    # [file path] specify output directory
                #"-A" : None,    # [file] import substitution alphabe
                "-L": 6,    # [integer] set minimum TR domain length
                "-P": 0,    # [0-1] set minimum % sequence coverage
                #"-Q" : None,    # [user,pass,db] send output to MySQL TR database
            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to detector"""
            toks = [
                optstring for optstring,
                optvalue in self.boolopts.items() if optvalue]

            toks = toks + \
                [optstring + str(optvalue)
                    for optstring, optvalue in self.valopts.items()
                    if optvalue is not None]

            if infile:
                toks.append(infile)

            return toks

    def __init__(self, name=name):
        """Construct DetectorXStream object.

        Arguments:
        name: The name of the tandem repeat detector.
        """

        self.config = self.Configuration()
        executable = BinaryExecutable(
            binary=REPEAT_DETECTOR_PATH[name],
            name=name)
        super(DetectorXStream, self).__init__(executable)

    def find_chartfile(self, searchdir):
        """Look for XSTREAM output xls file and return filepath"""
        files = os.listdir(searchdir)

        # look for the first file that matches XSTREAM_something_chart.xl
        pat_chartfilename = re.compile("^XSTREAM_.*_chart\.xls$")
        for f in files:
            if pat_chartfilename.match(f):
                return os.path.join(searchdir, f)
        else:
            return None

    def run_process(self, working_dir, infile):
        """ Run detector process on infile in working_dir and return all repeats found.

        Run detector process on infile in working_dir and return all repeats found.

        Args:
            working_dir (str): Working directory
            infile (str): Infile

        .. ToDo:: Complete docstrings.
        """

        wd = os.path.join(working_dir, DetectorXStream.name)

        # execute detector process
        retcode, stdoutfname, stderrfname = super(
            DetectorXStream, self).run_process(
            wd, *self.config.tokens(infile))

        if check_java_errors(stdoutfname, stderrfname,
                             log=LOG, procname=DetectorXStream.displayname):
            return []

        # Find output file, cry about it if not found

        #shutil.rmtree(os.path.join(EXECROOT, '..', 'spielwiese'))
        #shutil.copytree(working_dir, os.path.join(EXECROOT, '..', 'spielwiese'))

        chartfilename = self.find_chartfile(wd)
        if not chartfilename:
            LOG.error("No XSTREAM output chart file found in %s!", wd)
            return []

        # Process output file, return results
        with open(chartfilename, "r") as infile:
            return list(repeat_detection_io.xstream_get_repeats(infile))


def Detectors(lDetector=None, sequence_type=None):
    """ Define a global dictionary of all used detector functions.

    Define a global dictionary of all used detector functions.

    Args:
        lDetector (list of str): A list of repeat detection algorithm names.
        sequence_type (str): Either "AA" or "DNA".

    Raises:
        Exception: if at least one of the provided detectors in ``lDetector`` does not exist.

    """

    global DETECTORS

    if not sequence_type:
        sequence_type = CONFIG_GENERAL["sequence_type"]
    if not lDetector:
        lDetector = CONFIG_GENERAL["sequence"]["repeat_detection"][sequence_type]
    else:
        if isinstance(lDetector, str):
            lDetector = [lDetector]
        elif not isinstance(lDetector, list):
            raise TypeError(""" lDetector is not of type list. Please supply a list of TR detectors
            (e.g. lDetector = ['HHrepID']). If you use TR detectors defined in config.ini, make
            sure the TR detector list (AA or DNA) is correctly defined.""")
        if any(i not in FINDER_LIST.keys() for i in lDetector):
            raise Exception(
                "Unknown TR detector supplied (Supplied: {}. Known TR detectors: {})".format(
                    lDetector,
                    FINDER_LIST.keys()))

    DETECTORS = {FINDER_LIST[i].name: FINDER_LIST[i]() for i in lDetector}


def split_sequence(seq_records, working_dir):
    """ Split a FASTA file with multiple entries to several FASTA files with one entry

    Arguments:
    seq_records -- List of Sequence instances.
    working_dir -- Output directory for splitted file

    Returns:
    A list of tuples containing the Protein identifier and the file name

    .. TODO: Do not append i to outfiles, but a sequence ID
    """

    outfiles = []

    for i, seq in enumerate(seq_records):
        filename = "sequence_{0:03}.faa".format(i + 1)
        seq.write(file=os.path.join(working_dir, filename), file_format="fasta")
        outfiles.append((i, filename))

    return outfiles


class DetectorJob:

    def __init__(self, detector, infile, job_id, protein_id):
        self.detector = detector
        self.infile = infile
        self.job_id = job_id
        self.protein_id = protein_id


################ RUN A SET OF Tandem Repeat Detection algorithms (TRDs) ##

def run_detector(
        seq_records,
        detectors=None,
        sequence_type='AA',
        default=True,
        local_working_dir=None,
        num_threads=1):
    """ Run TRD on sequence_records and return the predicted repeats for each ``seq_records``
     and for each tandem repeat detector.

    Run TRD on sequence_records and return the predicted repeats for each ``seq_records``
    and for each tandem repeat detector. Example of return type::

        [
            # record 1
            {
            'T-REKS' : [ Repeat(), Repeat(), ...],
            'XSTREAM' : [ Repeat(), Repeat(), ...],
            ...
            },
            # record 2
            ...
        ]

    Args:
        seq_records (list of Sequence): A list of Sequence instances
        detectors (list of str): A list tandem repeat detector names
        sequence_type (str): Either "AA" or "DNA"
        default (bool): If True, default values for the detection algorithms are used.
        local_working_dir (str): Directory where data and results are stored. If provided,
        temporary files are not deleted/
        num_threads (int): Run ``num_threads`` detectors on parallel threads.

    Returns:
        list of dictionary: A list with a dictionary for each record in seq_records. The
        dictionary contains a list of repeats for each detector that was used.
    """

    # Create temporary working dir
    if local_working_dir:
        working_dir = local_working_dir
    else:
        working_dir = tempfile.mkdtemp()
        LOG.debug(
            "repeat_detection_run.run_detector: Created tempfile: %s",
            working_dir)
        if not os.path.isdir(working_dir):
            raise IOError("The specified directory \"" + working_dir +
                          "\" does not exist")

    # Initialise Detectors
    Detectors(detectors, sequence_type)

    # Adjust TRD parameters:
    if not default:
        if sequence_type == 'AA':
            DETECTORS['HHrepID'].config = set_hhrepid_config_open()
            DETECTORS['TRUST'].config = set_trust_config_open()
        else:
            DETECTORS['TRF'].config = set_trf_config_open()
            DETECTORS['PHOBOS'].config = set_phobos_config_open()
        DETECTORS['T-REKS'].config = set_treks_config_open()
        DETECTORS['XSTREAM'].config = set_xstream_config_open()

    if sequence_type == 'DNA':
        DETECTORS['T-REKS'].config = set_treks_config_DNA()

    infiles = split_sequence(seq_records, working_dir)

    # list of dictionaries for each protein containing detectorname : results
    # entries
    results = [
        {fname: [] for fname, detector in DETECTORS.items()}
        for i in range(len(infiles))
    ]

    # Create a list for our jobs
    joblist = [
        DetectorJob(
            detector=detector,
            infile=os.path.join(working_dir, infile[1]),
            job_id=job_id,
            protein_id=infile[0]
        )
        for fname, detector in DETECTORS.items()
        for job_id, infile in enumerate(infiles)
    ]

    LOG.info("Processing %d input files in %d jobs.",
             len(infiles), len(joblist))

    # put joblist into job queue, thus starting actual work
    for job in joblist:
        try:
            wd = os.path.join(working_dir, "{0:03}".format(job.job_id + 1))

            LOG.debug("Launching detector \"%s\" in directory %s",
                      job.detector.name, os.path.join(wd, job.detector.name))

            result = job.detector.run_process(wd, job.infile)

            # lock the mutex for results, append result
            # with result_lock:
            results[job.job_id][job.detector.name].extend(result)

            LOG.debug("Detector \"%s\" returned from job %d",
                      job.detector.name, job.job_id)
        except:
            LOG.exception(
                "Exception occured in worker while processing %s with %s",
                job.infile, job.detector.name)

    LOG.info("All jobs returned.")

    # delete temporary directory
    if not local_working_dir:
        try:
            shutil.rmtree(working_dir)
        # I guess the error type is known, and you could be more precise :)
        except:
            LOG.error("Unexpected error: {0}".format(sys.exc_info()[0]))
            raise

    return results


######## SET OPEN CONFIGS #########

def set_hhrepid_config_open():
    # construct open configuration for HHrepid
    config = DetectorHHrepID.Configuration()
    config.boolopts["-nofilt"] = True
    # <float> max p-value of suboptimal alignments in all search rounds but the last one (def=0.1)
    config.valopts["-P"] = 0.6
    config.valopts["-T"] = 0.2  # <float> max total repeat p-value (def=0.001)
    # [0,inf[  minimal length of repeats to be identified (def=7)
    config.valopts["-lmin"] = 0
    LOG.debug("%s config tokens: %s", DETECTORS['hhrepid'].displayname,
              ", ".join(DETECTORS['hhrepid'].config.tokens()))
    return config


def set_phobos_config_open():
    # construct open configuration for Phobos
    config = DetectorPhobos.Configuration()
    config.valopts["--maxUnitLen"] = 100
    config.valopts["--mismatchScore"] = -2
    config.valopts["--indelScore"] = -2
    return config

# Incomplete


def set_treks_config_DNA():
    # construct DNA configuration for T-Reks
    config = DetectorTREKS.Configuration()
    config.valopts["-type"] = 1
    LOG.debug("%s config tokens: %s", DETECTORS['t-reks'].displayname,
              ", ".join(DETECTORS['t-reks'].config.tokens()))
    return config


def set_treks_config_open():
    # construct open configuration for T-Reks
    config = DetectorTREKS.Configuration()
    config.valopts["-similarity"] = 0.2
    LOG.debug("%s config tokens: %s", DETECTORS['t-reks'].displayname,
              ", ".join(DETECTORS['t-reks'].config.tokens()))
    return config


def set_trf_config_open():
    # construct open configuration for TRF
    config = DetectorTRF.Configuration()
    config.valopts["Minscore"] = 20
    LOG.debug("%s config tokens: %s", DETECTORS['trf'].displayname,
              ", ".join(DETECTORS['trf'].config.tokens()))
    return config


def set_trust_config_open():
    # construct open configuration for Trust
    config = DetectorTrust.Configuration()
    config.boolopts['-force'] = True
    LOG.debug("%s config tokens: %s", DETECTORS['trust'].displayname,
              ", ".join(DETECTORS['trust'].config.tokens()))
    return config


def set_xstream_config_open():
    # construct open configuration for XStream
    config = DetectorXStream.Configuration()
    #config.boolopts[optname] = True
    config.valopts["-I"] = 0.2
    config.valopts["-i"] = 0.1
    config.valopts["-g"] = 100
    LOG.debug("%s config tokens: %s", DETECTORS['xstream'].displayname,
              ", ".join(DETECTORS['xstream'].config.tokens()))
    return config


######################## HARDCODED OVERVIEW DICTIONARIES #################

FINDER_LIST = {"HHrepID": DetectorHHrepID,
               "Phobos": DetectorPhobos,
               "TRED": DetectorTRED,
               "T-REKS": DetectorTREKS,
               "TRF": DetectorTRF,
               "TRUST": DetectorTrust,
               "XSTREAM": DetectorXStream
               }
