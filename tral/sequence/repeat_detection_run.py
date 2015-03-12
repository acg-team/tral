# (C) 2011, Alexander Korsunsky
# (C) 2011-2015 Elke Schaper

"""

    :synopsis: Execution of repeat detection algorithms

    .. moduleauthor:: Elke Schaper <elke@inf.ethz.ch>

"""

import itertools
import logging
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile

from collections import OrderedDict

from tral import configuration
from tral.sequence import repeat_detection_io
from tral.paths import *

log = logging.getLogger(__name__)

c = configuration.Configuration.Instance()
general_config = c.config
repeat_detector_path = general_config["sequence"]["repeat_detector_path"]

MAX_MEMORY_USAGE = str(10000000)

class BinaryExecutable:
    def __init__(self, binary=None):
        """Construct a BinaryExecutable object.

        """

        if not binary:
            raise TypeError("A binary executable must be provided :) ")
        self.binary = binary

    def get_execute_tokens(self, *args):
        """Return the tokens to invoke the program with the arguments args"""

        return [self.binary] + list(args)


    def get_execute_line(self, *args):
        """Return the command line to invoke the program with the arguments args"""
        return " ".join(self.get_execute_tokens(*args))


class TRFFinder(object):

    def __init__(self, executable):
        """Construct a TRFFinder object with executable als executable object"""
        log.debug(executable)
        self.__executable = executable

    def run_process(self, working_dir, *args):
        """Launch a finder process

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

        log.debug("Launching process: %s in %s",
            self.__executable.get_execute_line(*args), working_dir)

        log.debug("Launching process tokens: %s in %s",
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


class FinderHHrepID(TRFFinder):
    name = 'HHrepID'
    displayname = "HHrepID"

    """ Execute vi
    ./hhrepid_32 -i infile -v 0 -nofilt -d dummyHMM.hmm -o resultfile
    ./hhrepid -i <query> -d <path to cal.hhm> -tp <path to tp.dat> -fp <path to fp.dat> """

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "-nofilt"        : False,    # turn off low-complexity filter (default: filter on)
                "-shuffle"       : False    # calibration by shuffling instead of database search (default: off)
            }

            self.valopts = {
                "-i"      :None,        # <file> input query alignment  (fasta/a2m/a3m) or HMM file (.hhm)
                "-d"      :repeat_detector_path['HHrepID_dummyhmm'],   # <path> dummy hmm database file
                "-o"      :'hhrepID.o',    # <file> write results and multiple sequence alignment to file (default=none)
                "-v"      :0,           # -v: verbose mode (default: show only warnings)  ;  -v 0: suppress all screen outpu
                "-P"      :None,        # <float> max p-value of suboptimal alignments in all search rounds but the last one (def=0.1)
                "-R"      :None,        # <float> max p-value of repeats  (def=1)
                "-T"      :None,        # <float> max total repeat p-value (def=0.001)
                "-alpha"  :None,        # <float> For calculating repeat p-values: weight of n'th suboptimal alignment vs. (n+1)-st suboptimal alignmen
                "-k"      :None,        # [0,inf[ For calculating repeat p-values: maximal number of suboptimal alignments considered (def=1)
                "-cont"   :None,        # <float>  probability threshold for masking of inconsistent cells (def=0.001)
                "-mrgr"   :None,        # [0,inf[  number of merge rounds for achieving consistency (def=3)
                "-lcon"   :None,        # <float>  preserved local context in shuffling (def=2)
                "-lmin"   :None,        # [0,inf[  minimal length of repeats to be identified (def=7)

            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""
            toks = []
            if infile:
                toks += ['-i', infile]

            for optstring, optvalue in self.valopts.items():
                if optvalue != None:
                    toks += [str(optstring), str(optvalue)]

            for optstring, optvalue in self.boolopts.items():
                if optvalue: toks.append(optstring)

            return toks


    def __init__(self,
        executable=BinaryExecutable(binary = repeat_detector_path[name]),
        config = Configuration()
    ):
        """Construct FinderHHrepID object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderHHrepID, self).__init__(executable)
        if config == None:
            self.config = FinderHHrepID.Configuration()
        else:
            self.config = config

    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, FinderHHrepID.name)

        prog_args = self.config.tokens(infile=infile)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args) ##

        ## CHECK WHETHER outfile exists, else give warning and return empty list.

        # Process output file, return results
        resultfile = os.path.join(wd, self.config.valopts["-o"])
        if os.path.isfile(resultfile):
            with open(resultfile, "r") as infilehandle:
                tmp = list(repeat_detection_io.hhpredid_get_repeats(infilehandle))
            return tmp
        else:
            return []



class FinderPhobos(TRFFinder):
    name = 'PHOBOS'
    displayname = "PHOBOS"

    # execute via phobos --printRepeatSeqMode 3 DNA.faa phobos.o
    # or phobos --help to see all options

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "--preferShorterRepeats" : False,   # Phobos will choose the shorter of the two alignments instead of the longer. Phobos will run faster with this option.
                "--dontRemoveMostlyOverlapping": False, # dontRemoveMostlyOverlapping. Default: DO remove mostly overlapping repeats in favour of the one with the higher internal score
            }

            self.scriptopts = {
                "outputfile"            : 'phobos.o',
            }

            self.valopts = {
                "--reportUnit" :  None, # <int> Display of TR. 0: asIs, 1: Alphabetical normal form, 2: Alphabetical normal form also considering the reverse complement. Default: 2
                "--NPerfectionMode" : None, # <int> Treatment of N.  0: asMismatch, 1: asNeutral, 2: asMatch. Default: 0.
                "--printRepeatSeqMode": '3', # <int> Show the repeat unit alignment
                "--outputFormat": None, # <int> 0 : Use the Phobos output format. Default: 0
                "--minPerfection": None, # <float> Minimum perfection of a satellites. Default: 0.
                "--maxPerfection": None, # <float> Maximum perfection of a satellites. Default: 100.
                "--maximum_score_reduction": None, # The maximum amount the score can be reduced before search is aborded. Typical: 6*mismatch-penalty or infinite. Default: infinite
                "--recursion": None, # The recursion depth used in the search. Values in the range 3 to 7 are recommended. A value of 0 implies a search for perfect repeats only. Default: Not stated
                "--minUnitLen": None, # <int> Minimum unit length. Default: 1
                "--maxUnitLen": None, # <int> Maximum unit length. Default: 10
                "--minLength_b": None, # <float> The minimum length of a repeat is determined with: maximum( minLength, minLength_a + minLength_b*(unit-length) ). Default value of minLength_b: 0
                "--minLength_a": None, # <int> Default value of minLength_a: 0
                "--minLength": None, # <int> Default value of minLength: 0
                "--minScore_b": None, # <float> The minimum score of a repeat is determined with: maximum( minScore, minScore_a + minScore_b*(unit-length) ). Default value of minScore_b: 1
                "--minScore_a": None, # <int> Default value of minScore_a: 0
                "--minScore": None, # <int> Default value of minScore: 6
                "--mismatchScore": None, # <int> Score for mismatch - must be negative. Default: -5. Match score is fixed to one.
                "--indelScore": None, # <int> Score for indels - must be negative. Default: -5. Match score is fixed to one.
                "--searchMode": 'imperfect' # ['exact','extendExact','imperfect']
            }


        def set_working_dir(self, working_dir):
            self.scriptopts['outputfile'] = os.path.join(working_dir,self.scriptopts['outputfile'])
            if not os.path.isdir(working_dir):
                os.makedirs(working_dir)

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""
            toks = [optstring
                for optstring, optvalue in self.boolopts.items() if optvalue]

            for optstring, optvalue in self.valopts.items():
                if optvalue != None:
                    toks.append(optstring)
                    toks.append(str(optvalue))

            if infile:
                toks.append(infile)
            toks.append(self.scriptopts['outputfile'])
            return toks

    def __init__(self,
        executable=BinaryExecutable(binary=os.path.join(repeat_detector_path[name],"phobos")),
        config = None
    ):
        """Construct FinderPhobos object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """

        super(FinderPhobos, self).__init__(executable)
        if config == None:
            self.config = FinderPhobos.Configuration()
        else:
            self.config = config

    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, FinderPhobos.name)
        self.config.set_working_dir(working_dir=wd)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
               super(FinderPhobos, self).run_process(wd, *self.config.tokens(infile))
        ## alternatively:
        #prog_args = self.config.tokens(infile=infile)
        #retcode, stdoutfname, stderrfname = super().run_process(wd, *prog_args)

        # Process output file, return results If outfile does not exist return empty list
        if os.path.isfile(self.config.scriptopts['outputfile']):
            with open(self.config.scriptopts['outputfile'], "r") as resultfilehandle:
                tmp = list(repeat_detection_io.phobos_get_repeats(resultfilehandle))
            return tmp
        else:
            log.warning("Did not find Phobos result file in %s", self.config.scriptopts['outputfile'])
            return []


class FinderTRED(TRFFinder):
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

            infile -- generate tokens to pass infile as input file to finder"""
            toks = []
            if infile:
                toks = [infile]
            toks.append(self.result_file)

            return toks

    def __init__(self,
        executable=BinaryExecutable(binary=repeat_detector_path[name]),
        config = None
    ):
        """Construct FinderTRED object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderTRED, self).__init__(executable)
        if config == None:
            self.config = FinderTRED.Configuration()
        else:
            self.config = config

    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, FinderTRED.name)

        self.result_file = os.path.join(wd, 'tred.o')

        self.config.set_result_file(result_file=self.result_file)
        prog_args = self.config.tokens(infile=infile)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args) ##

        ## CHECK WHETHER outfile exists, else give warning and return empty list.

        # Process output file, return results
        if os.path.isfile(self.result_file):
            with open(self.result_file, "r") as infilehandle:
                tmp = list(repeat_detection_io.tred_get_repeats(infilehandle))
            return tmp
        else:
            log.warning("Did not find TRED result file in %s", self.result_file)
            return []


class FinderTREKS(TRFFinder):
    name = 'T-REKS'
    displayname = "T-REKS"

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "-overlapfilter" : True,
                "-nosplit" : True
            }

            self.valopts = {
                "-msaMode"      :None,
                "-clustal"      :None,
                "-db"           :None,
                "-muscle"       :None,
                "-outfile"      :None,
                "-similarity"   :None,
                "-kmeans"       :None,
                "-varIndels"    :None,
                "-type"         :None
            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""
            bool_toks = [optstring
                for optstring, optvalue in self.boolopts.items() if optvalue]
            value_toks = ([
                optstring+"="+str(optvalue)
                for optstring, optvalue in self.valopts.items()
                if optvalue != None
            ])

            if infile:
                value_toks.append("-infile="+infile)

            return bool_toks + value_toks

    def __init__(self,
        executable=BinaryExecutable(binary = repeat_detector_path[name]),
        config = None
    ):
        """Construct FinderTREKS object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderTREKS, self).__init__(executable)
        if config == None:
            self.config = FinderTREKS.Configuration()
        else:
            self.config = config


    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, FinderTREKS.name)

        prog_args = self.config.tokens(infile=infile)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super(FinderTREKS, self).run_process(wd, *prog_args)

        if check_java_errors(stdoutfname, stderrfname,
            log=log, procname=FinderTREKS.displayname):
            return []

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            tmp = list(repeat_detection_io.treks_get_repeats(outfile))
        return tmp

class FinderTRF(TRFFinder):
    name = 'TRF'
    displayname = "TRF_Benson"

    ### Execute via:
    ## ./trf404.linux64.exe inputfile 2 7 7 80 10 50 500 -d
    ## Result file in this case is called inputfile.2.7.7.80.10.50.500.1.txt.html
    ## WARNING: This result file name is only used if there is only one sequence in the input file!
    ## result_file_name = ".".join([inputfile] + [str(optvalue) for optstring, optvalue in self.valopts.items() if optvalue != None] + ['1.txt.html'])

    ## ./trf404.linux64.exe File Match Mismatch Delta PM PI Minscore MaxPeriod [options]

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "-r"            : False,    # no redundancy elimination # Already without the -r flag set, the three best redundant solutions are displayed
                #"-h"            : False,    # When -h is set, the output html of interest is not written :(
                #"-d"            : False,    # produce data file, not needed at curren
                #"-m"            : False,    # masked sequence file
                #"-f"            : False,    # flanking sequence
            }

            ## The order of TRF flags matters. Hence, we use an ordered dictionary for valopts.
            valopts_tmp = [
                ("File"          ,None),  # input file
                ("Match"         ,2),     # matching weigh
                ("Mismatch"      ,7),     # mismatching penalty
                ("Delta"         ,7),     # indel penalty
                ("PM"            ,80),    # integer match probability
                ("PI"            ,10),    # integer indel probability
                ("Minscore"      ,50),    # minimum alignment score to report
                ("MaxPeriod"     ,1000)    # maximum period size to report

            ]
            self.valopts = OrderedDict(valopts_tmp)

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""

            toks = []
            if infile:
                toks.append(infile)

            for optstring, optvalue in self.valopts.items():
                if optvalue != None:
                    toks.append(str(optvalue))

            for optstring, optvalue in self.boolopts.items():
                if optvalue: toks.append(optstring)

            return toks


    def __init__(self,
        executable=BinaryExecutable(binary=repeat_detector_path[name]),
        config = None
    ):
        """Construct FinderTRF object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderTRF, self).__init__(executable)
        if config == None:
            self.config = FinderTRF.Configuration()
        else:
            self.config = config


    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""

        wd = os.path.join(working_dir, FinderTRF.name)

        prog_args = self.config.tokens(infile=infile)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super().run_process(wd, *prog_args) ##

        # get the infile name from the infile path
        infile_name = re.findall(r"/([_\.\w]+)$", infile)[0]

        ### WARNING! TRF is a volatile programme, this result_file_name might not always be true.
        ### Maybe a general search on txt.html elements in the result dir might be more appropiate.
        result_file_name = ".".join([infile_name] + [str(optvalue) for optstring, optvalue in self.config.valopts.items() if optvalue != None] + ['1.txt.html'])

        # Process output file, return results
        with open(os.path.join(wd, result_file_name), "r") as outfile:
            tmp = list(repeat_detection_io.trf_get_repeats(outfile))
        return tmp



class FinderTrust(TRFFinder):
    name = 'TRUST'
    displayname = "TRUST"

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "-noseg" : True,
                "-force" : False,
            }

            self.valopts = {
                "-matrix" : repeat_detector_path['TRUST_substitutionmatrix'],
                "-gapo" : "8",
                "-gapx" : "2",
                "-procTotal": "1",
            } # procTotal: When running on a cluster, specify the amount of processors used  # procNr cpu_number:  When running on a cluster, specify 0-based CPU-number and amount of processors used

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""
            toks = [optstring
                for optstring, optvalue in self.boolopts.items() if optvalue]

            for optstring, optvalue in self.valopts.items():
                if optvalue != None:
                    toks.append(optstring)
                    toks.append(str(optvalue))

            if infile:
                toks.append("-fasta")
                toks.append(infile)

            return toks

    def __init__(self,
        executable=BinaryExecutable(binary = repeat_detector_path[name]),
        config = Configuration()
    ):
        """Construct FinderTrust object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderTrust, self).__init__(executable)
        self.config = config


    def run_process(self, working_dir, infile):
        """Run finder process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, FinderTrust.name)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
               super(FinderTrust, self).run_process(wd, *self.config.tokens(infile))

        if check_java_errors(stdoutfname, stderrfname,
            log=log, procname=FinderTrust.displayname):
            return []

        #shutil.copyfile(stdoutfname, YOUR FAVOURITE PATH)

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            return list(repeat_detection_io.trust_get_repeats(outfile))


class FinderXStream(TRFFinder):
    name = 'XSTREAM'
    displayname = "XSTREAM"

    class Configuration:
        def __init__(self):
            self.boolopts = {
                #"-n" : False,   # sequence type: nucleotide
                "-l" : False,   # don't look for periods >1000
                #"-p" : False,   # don't print repeat alignmen
                "-h" : True,    # write to console, not HTML
                #"-C" : False,   # don't print any repeat info
                "-z" : True,  # generate spreadsheet .xls output
                #"-Z" : False,   # create statistics spreadshee
                #"-s" : False,   # parse output by input sequence
                "-f" : False,   # perform divide & conquer on inpu
                #"-G" : False,   # create color-coded TR block diagram
                #"-G*" : False,  # create single color TR block diagram
                #"-B" : False,   # create sequence comparison map PNG
                #"-O" : False,   # print "mismatch/gaps colored" HTML outpu
                "-N" : False,   # turn off nesting
                "-o" : False,   # don't remove overlapping TR domain
                #"-S" : False,   # input strain into database
                "-V" : True,   # don't invoke D & C automatically for long inpu
            }

            self.valopts = {
                "-g" : None,    # [integer] set gap number
                "-i" : None,    # [0-1] set matching threshold for detection
                "-I" : None,    # [0-1] set cons error threshold for printing
                "-D" : None,    # [0-1] set max % indels for printing
                "-m" : 1,    # [integer] set minimum period
                "-x" : None,    # [integer] set maximum period
                "-e" : None,    # [any number>=2] set minimum copy number
                # "-a" : None,    # [string] add identifier string to output file
                "-f" : None,    # [integer] set fragment length
                #"-d" : None,    # [file path] specify output directory
                #"-A" : None,    # [file] import substitution alphabe
                "-L" : 6,    # [integer] set minimum TR domain length
                "-P" : 0,    # [0-1] set minimum % sequence coverage
                #"-Q" : None,    # [user,pass,db] send output to MySQL TR database
            }

        def tokens(self, infile=None):
            """Generate command line tokens based on this configuration.
            Arguments:

            infile -- generate tokens to pass infile as input file to finder"""
            toks = [optstring
                for optstring, optvalue in self.boolopts.items() if optvalue]

            toks = toks + \
                [optstring + str(optvalue)
                    for optstring, optvalue in self.valopts.items()
                    if optvalue != None]

            if infile:
                toks.append(infile)

            return toks

    def __init__(self,
        executable=BinaryExecutable(binary = repeat_detector_path[name]),
        config = Configuration()
    ):
        """Construct FinderXStream object.

        Arguments:
        executable -- Use this executable object instead of default-constructed one
        config -- Use this configuration object instead of default-constructed one
        """
        super(FinderXStream, self).__init__(executable)
        self.config = config

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
        """ Run finder process on infile in working_dir and return all repeats found.

        Run finder process on infile in working_dir and return all repeats found.

        Args:
            working_dir (str): Working directory
            infile (str): Infile

        .. ToDo:: Complete docstrings.
        """

        wd = os.path.join(working_dir, FinderXStream.name)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super(FinderXStream, self).run_process(wd, *self.config.tokens(infile))


        if check_java_errors(stdoutfname, stderrfname,
            log=log, procname=FinderXStream.displayname):
            return []

        # Find output file, cry about it if not found

        #shutil.rmtree(os.path.join(EXECROOT, '..', 'spielwiese'))
        #shutil.copytree(working_dir, os.path.join(EXECROOT, '..', 'spielwiese'))

        chartfilename = self.find_chartfile(wd)
        if not chartfilename:
            log.error("No XSTREAM output chart file found in %s!", wd)
            return []

        # Process output file, return results
        with open(chartfilename, "r") as infile:
            return list(repeat_detection_io.xstream_get_repeats(infile))



def Finders(lFinder = None, sequence_type = None):
    """ Define a global dictionary of all used finder functions.

    Define a global dictionary of all used finder functions.

    Args:
        lFinder (list of str): A list of repeat detection algorithm names.
        sequence_type (str): Either "AA" or "DNA".

    Raises:
        Exception: if at least one of the provided finders in ``lFinder`` does not exist.

    .. ToDo:: Is FINDER_LIST defined correctly?
    """

    global finders

    if not sequence_type:
        sequence_type = general_config["sequence_type"]
    if not lFinder:
        lFinder = general_config["sequence"]["repeat_detection"][sequence_type]
    else:
        if any(i not in list(itertools.chain(*general_config["sequence"]["repeat_detection"].values())) for i in lFinder):
            raise Exception("Unknown TR detector supplied (Supplied: {}. Known TR detectors: {})".format(lFinder, FINDER_LIST))

    finders = {FINDER_LIST[i]:FINDER_FUNCTION_LIST[i] for i in lFinder}


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
        filename = "sequence_{0:03}.faa".format(i+1)
        seq.write(file = os.path.join(working_dir,filename), format = "fasta")
        outfiles.append((i, filename))

    return outfiles


class FinderJob:
    def __init__(self, finder, infile, job_id, protein_id):
        self.finder = finder
        self.infile = infile
        self.job_id = job_id
        self.protein_id = protein_id



################ RUN A SET OF Tandem Repeat Detection algorithms (TRDs) ##################

def run_TRD(seq_records, lFinders = None, sequence_type = 'AA', default = True, local_working_dir=None, num_threads=1):
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
        lFinders (list of str): A list tandem repeat detector names
        sequence_type (str): Either "AA" or "DNA"
        default (bool): If True, default values for the detection algorithms are used.
        local_working_dir (str): Directory where data and results are stored. If provided,
        temporary files are not deleted/
        num_threads (int): Run ``num_threads`` finders on parallel threads.

    Returns:
        list of dictionary: A list with a dictionary for each record in seq_records. The
        dictionary contains a list of repeats for each finder that was used.
    """

    # Create temporary working dir
    if local_working_dir:
        working_dir = local_working_dir
    else:
        working_dir = tempfile.mkdtemp()
        log.debug("repeat_detection_run.run_TRD: Created tempfile: %s", working_dir)
        if not os.path.isdir(working_dir):
            raise IOError("The specified directory \""+working_dir+
                          "\" does not exist")

    # Initialise Finders
    Finders(lFinders, sequence_type)

    ## Adjust TRD parameters:
    if not default:
        if sequence_type == 'AA':
            finders['HHrepID'].config = set_hhrepid_config_open()
            finders['TRUST'].config = set_trust_config_open()
        else:
            finders['TRF'].config = set_trf_config_open()
            finders['PHOBOS'].config = set_phobos_config_open()
        finders['T-REKS'].config = set_treks_config_open()
        finders['XSTREAM'].config = set_xstream_config_open()

    if sequence_type == 'DNA':
        finders['T-REKS'].config = set_treks_config_DNA()


    infiles = split_sequence(seq_records, working_dir)

    # list of dictionaries for each protein containing findername : results
    # entries
    results = [
        {fname : [] for fname,finder in finders.items()}
            for i in range(len(infiles))
    ]

    # Create a list for our jobs
    joblist = [
        FinderJob(
            finder=finder,
            infile=os.path.join(working_dir,infile[1]),
            job_id=job_id,
            protein_id=infile[0]
        )
        for fname, finder in finders.items()
        for job_id, infile in enumerate(infiles)
    ]

    log.info("Processing %d input files in %d jobs.",
        len(infiles), len(joblist))

    # put joblist into job queue, thus starting actual work
    for job in joblist:
        try:
            wd = os.path.join(working_dir, "{0:03}".format(job.job_id+1))

            log.debug("Launching finder \"%s\" in directory %s",
                job.finder.name, os.path.join(wd, job.finder.name))

            result = job.finder.run_process(wd, job.infile)

            # lock the mutex for results, append result
            #with result_lock:
            results[job.job_id][job.finder.name].extend(result)

            log.debug("Finder \"%s\" returned from job %d",
                job.finder.name, job.job_id)
        except:
            log.exception(
                "Exception occured in worker while processing %s with %s",
                job.infile, job.finder.name)

    log.info("All jobs returned.")

    # delete temporary directory
    if not local_working_dir:
        try:
            shutil.rmtree(working_dir)
        except: # I guess the error type is known, and you could be more precise :)
            logging.error("Unexpected error: {0}".format(sys.exc_info()[0]))
            raise

    return results


######## SET OPEN CONFIGS #########

def set_hhrepid_config_open():
    # construct open configuration for HHrepid
    config = FinderHHrepID.Configuration()
    config.boolopts["-nofilt"] = True
    config.valopts["-P"] = 0.6  # <float> max p-value of suboptimal alignments in all search rounds but the last one (def=0.1)
    config.valopts["-T"] = 0.2  # <float> max total repeat p-value (def=0.001)
    config.valopts["-lmin"] = 0  # [0,inf[  minimal length of repeats to be identified (def=7)
    log.debug("%s config tokens: %s", finders['hhrepid'].displayname,
                        ", ".join(finders['hhrepid'].config.tokens()))
    return config

def set_phobos_config_open():
    # construct open configuration for Phobos
    config = FinderPhobos.Configuration()
    config.valopts["--maxUnitLen"] = 100
    config.valopts["--mismatchScore"] = -2
    config.valopts["--indelScore"] = -2
    return config

## Incomplete
def set_treks_config_DNA():
    # construct DNA configuration for T-Reks
    config = FinderTREKS.Configuration()
    config.valopts["-type"] = 1
    log.debug("%s config tokens: %s", finders['t-reks'].displayname,
                        ", ".join(finders['t-reks'].config.tokens()))
    return config

def set_treks_config_open():
    # construct open configuration for T-Reks
    config = FinderTREKS.Configuration()
    config.valopts["-similarity"] = 0.2
    log.debug("%s config tokens: %s", finders['t-reks'].displayname,
                        ", ".join(finders['t-reks'].config.tokens()))
    return config

def set_trf_config_open():
    # construct open configuration for TRF
    config = FinderTRF.Configuration()
    config.valopts["Minscore"] = 20
    log.debug("%s config tokens: %s", finders['trf'].displayname,
                        ", ".join(finders['trf'].config.tokens()))
    return config

def set_trust_config_open():
    # construct open configuration for Trust
    config = FinderTrust.Configuration()
    config.boolopts['-force'] = True
    log.debug("%s config tokens: %s", finders['trust'].displayname,
                        ", ".join(finders['trust'].config.tokens()))
    return config


def set_xstream_config_open():
    # construct open configuration for XStream
    config = FinderXStream.Configuration()
    #config.boolopts[optname] = True
    config.valopts["-I"] = 0.2
    config.valopts["-i"] = 0.1
    config.valopts["-g"] = 100
    log.debug("%s config tokens: %s", finders['xstream'].displayname,
                        ", ".join(finders['xstream'].config.tokens()))
    return config


######################## HARDCODED OVERVIEW DICTIONARIES #################################

FINDER_FUNCTION_LIST = { "HHrepID": FinderHHrepID(),
                "Phobos": FinderPhobos(),
                "TRED": FinderTRED(),
                "T-REKS": FinderTREKS(),
                "TRF": FinderTRF(),
                "TRUST": FinderTrust(),
                "XSTREAM": FinderXStream()
                }

FINDER_LIST = { "HHrepID": FinderHHrepID.name,
                "Phobos": FinderPhobos.name,
                "TRED": FinderTRED.name,
                "T-REKS": FinderTREKS.name,
                "TRF": FinderTRF.name,
                "TRUST": FinderTrust.name,
                "XSTREAM": FinderXStream.name
                }