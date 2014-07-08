# (C) 2011, Alexander Korsunsky
# (C) 2011-2014 Elke Schaper

import itertools
import logging
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading

from collections import OrderedDict
from Bio import SeqIO

from tandemrepeats.sequence import repeat_detection_io
from tandemrepeats.paths import *

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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


class JavaExecutable:
    def __init__(self, javaclass=None, classpaths=[], jars=[],
                 javabin="java", javaopts=[]):
        """Construct a JavaExecutable object.

        Arguments:
        javaclass -- A string denoting the java class to be loaded for the
            program entry point. If None, don't load any java classe
        classpaths -- List of directories to search for java classe
        jars -- List of jars to search for java classe
        javabin -- java executable location
        javaopts -- Additional java option

        Either javaclass or at least one jar have to be provided.
        """

        if not javaclass and len(jars)==0:
            raise TypeError("Either a javaclass or a jar file\
                must be provided!")

        # assign object attribute
        self.__javabin = javabin
        self.__javaclass = javaclass

        # prepend "-jar" or "-cp" to all jars and classpath
        jaropts = ["-jar" for i in range(2*len(jars)) ]
        jaropts[1:len(jaropts):2] = jars
        cp_opts = ["-cp" for i in range(2*len(classpaths))]
        cp_opts[1:len(cp_opts):2] = classpaths

        self.__javaopts = javaopts + cp_opts + jaropts


    def get_execute_tokens(self, *args):
        """Return the tokens to invoke the program with the arguments args"""

        # only return the Java class if it is not empty
        if self.__javaclass:
            classes = [self.__javaclass]
        else:
            classes = []

        return [self.__javabin] + self.__javaopts + classes + list(args)


    def get_execute_line(self, *args):
        """Return the command line to invoke the program with the arguments args"""
        return " ".join(self.get_execute_tokens(*args))



def check_java_errors(outfile, errfile, logger=None, procname=None):
    """Checks for java problems and returns True if there were problem
    and False if there weren't.

    Arguments:
    outfile  -- Name of the redirected standard output channel file
    errfile  -- Name of the redirected standard error channel file
    logger   -- Name of the logger to issue log messages to. If none, no log
                messages will be issued.


    Checks for the following things:
     - Stdout file is empty but stderr file is not.
     - Java Exception string is indicated in the errfile
    """

    has_error = False

    with open(outfile, mode='r') as of, open(errfile, mode='r') as ef:
        outfile_size = of.seek(0, os.SEEK_END)
        errfile_str = ef.read()

    errfile_size = len(errfile_str)

    if (errfile_size != 0 and outfile_size == 0):
        has_error = True
        if logger:
            logger.warning(
                "Process \"%s\" has empty STDOUT but non-empty STDERR",
                procname)

    pat_javaexc = re.compile(r'Exception in thread ".+" \S+:.*$(\s+at .+)+', re.M)

    m = pat_javaexc.search(errfile_str)
    if m:
        has_error = True
        if logger:
            logger.warning(
                "Java Exception probably occured in process \"%s\"! Exception information:\n%s",
                procname, m.group(0))

    return has_error


class TRFFinder(object):

    def __init__(self, executable):
        """Construct a TRFFinder object with executable als executable object"""
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
        be redirected to working_dir/stdout.txt and working_dir/stdout.tx
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

        logger.debug("Launching process: %s in %s",
            self.__executable.get_execute_line(*args), working_dir)

        logger.debug("Launching process tokens: %s in %s",
            self.__executable.get_execute_tokens(*args), working_dir)
        # launch process
        __process = subprocess.Popen(
            self.__executable.get_execute_tokens(*args),
            cwd=working_dir, stdout=__stdout_file, stderr=__stderr_file, close_fds=True
        )

        __process.wait()

        __stdout_file.close()
        __stderr_file.close()

        return __process.returncode, stdoutfname, stderrfname


class FinderHHrepID(TRFFinder):
    name = 'hhrepid' # or 'hhrepID_64' ?
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
                "-d"      :os.path.join(DATAROOT,'hhrepid','dummyHMM.hmm'),   # <path> dummy hmm database file
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
        executable=BinaryExecutable(binary=os.path.join(EXECROOT,"hhrepid_64")),
        config = None
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
    name = 'phobos'
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
        executable=BinaryExecutable(binary=os.path.join(EXECROOT,"phobos")),
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
            logger.warning("Did not find Phobos result file in %s", self.config.scriptopts['outputfile'])
            return []


class FinderTRED(TRFFinder):
    name = 'tred'
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
        executable=BinaryExecutable(binary=os.path.join(EXECROOT,"tred")),
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
            logger.warning("Did not find TRED result file in %s", self.result_file)
            return []


class FinderTREKS(TRFFinder):
    name = 't-reks'
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
        executable=JavaExecutable(javaopts = ['-mx512m'], jars=[os.path.join(EXECROOT,"T-Reks.jar")]),
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
            logger=logger, procname=FinderTREKS.displayname):
            return []

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            tmp = list(repeat_detection_io.treks_get_repeats(outfile))
        return tmp

class FinderTRF(TRFFinder):
    name = 'trf'
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
        executable=BinaryExecutable(binary=os.path.join(EXECROOT,"trf")),
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
    name = 'trust'
    displayname = "TRUST"

    class Configuration:
        def __init__(self):
            self.boolopts = {
                "-noseg" : True,
                "-force" : False,
            }

            self.valopts = {
                "-matrix" : os.path.join(EXECROOT, "trust", "BLOSUM50"),
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
        executable=JavaExecutable(javaopts = ['-mx512m'],
            javaclass='nl.vu.cs.align.SelfSimilarity',
            classpaths=[os.path.join(EXECROOT,"trust")]
        ),
        config=Configuration()
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
            logger=logger, procname=FinderTrust.displayname):
            return []

        #shutil.copyfile(stdoutfname, YOUR FAVOURITE PATH)

        # Process output file, return results
        with open(stdoutfname, "r") as outfile:
            return list(repeat_detection_io.trust_get_repeats(outfile))


class FinderXStream(TRFFinder):
    name = 'xstream'
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
        executable=JavaExecutable(javaopts = ['-mx1024m'], jars=[os.path.join(EXECROOT,"xstream.jar")]),
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
        """Run finder process on infile in working_dir and return all repeats found"""
        wd = os.path.join(working_dir, FinderXStream.name)

        # execute finder process
        retcode, stdoutfname, stderrfname = \
            super(FinderXStream, self).run_process(wd, *self.config.tokens(infile))


        if check_java_errors(stdoutfname, stderrfname,
            logger=logger, procname=FinderXStream.displayname):
            return []

        # Find output file, cry about it if not found

        #shutil.rmtree(os.path.join(EXECROOT, '..', 'spielwiese'))
        #shutil.copytree(working_dir, os.path.join(EXECROOT, '..', 'spielwiese'))

        chartfilename = self.find_chartfile(wd)
        if not chartfilename:
            logger.error("No XSTREAM output chart file found in %s!", wd)
            return []

        # Process output file, return results
        with open(chartfilename, "r") as infile:
            return list(repeat_detection_io.xstream_get_repeats(infile))



def Finders(lFinder = None, sequence_type = "AA"):

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

    if not lFinder:
        lFinder = FINDER_DEFAULT[sequence_type]
    else:
        if any(i not in list(itertools.chain(*FINDER_DEFAULT.values())) for i in lFinder):
            raise Error("Unknown TR detector supplied (Supplied: {}. Known TR detectors: {})".format(lFinder, FINDERLIST))

    finders = {FINDER_LIST[i]:FINDER_FUNCTION_LIST[i] for i in lFinder}


def split_sequence(seq_records, working_dir):
    """ Split a FASTA file with multiple entries to several FASTA files with one entry

    Arguments:
    seq_records -- Iterator to the Bio.SeqIO record objec
    working_dir -- Output directory for splitted file

    Returns:
    A list of tuples containing the Protein identifier and the file name
    """
    # FIXME return protein length here!

    outfiles = []


    for i, record in enumerate(seq_records):
        filename = "sequence_{0:03}.faa".format(i+1)
        with open(os.path.join(working_dir,filename), "w", encoding="UTF-8") as outhandle:
            SeqIO.write(record, outhandle, "fasta")
        outfiles.append((record.id, filename))

    return outfiles


class FinderJob:
    def __init__(self, finder, infile, job_id, protein_id):
        self.finder = finder
        self.infile = infile
        self.job_id = job_id
        self.protein_id = protein_id

def finder_worker(working_dir, job_queue, result_list, result_lock):
    def clear_queue(q):
        """Eat all items in job queue.
        To be called if you wish to terminate execution."""
        while not q.empty():
            try:
                q.get_nowait()
                q.task_done()
            except:
                break

    # Run infinitely, because this should be called as a daemon thread
    while True:
        job = job_queue.get()

        try:
            wd = os.path.join(working_dir, "{0:03}".format(job.job_id+1))

            logger.debug("Launching finder \"%s\" in directory %s",
                job.finder.name, os.path.join(wd, job.finder.name))

            result = job.finder.run_process(wd, job.infile)


            # lock the mutex for results, append resul
            with result_lock:
                result_list[job.job_id][job.finder.name].extend(result)


            logger.debug("Finder \"%s\" returned from job %d",
                job.finder.name, job.job_id)
        except:
            # On errors, kill all jobs in queue, cry to the logfile and return
            clear_queue(job_queue)
            logger.exception(
                "Exception occured in worker while processing %s with %s",
                job.infile, job.finder.name)
            # raise
        finally:
            job_queue.task_done()


def run_finders(seq_records, working_dir=None, num_threads=1):
    """ Run all finder modules.

    Keyword arguments:

    seq_records: An iterator to Bio.SeqIO sequence record
    working_dir: operate in this directory.
        If None, creates a new temporary directory.
    num_threads: run num_threads finders simultanously

    Returns:

    A list with a dictionaryfor each record in seq_records. The dictionary
    contains a list of repeats for each finder that was used.
    example:
    [
        # record 1
        {
        't-reks' : [ Repeat(), Repeat(), ...],
        'xstream' : [ Repeat(), Repeat(), ...],
        ...
        },
        # record 2
        ...
    ]
    """

    # create temporary directory if none was specified
    if working_dir == None:
        working_dir = tempfile.mkdtemp()
    else:
        working_dir = os.path.abspath(working_dir)

    if not os.path.isdir(working_dir):
        raise IOError("The specified directory \""+working_dir+
                      "\" does not exist")

    infiles = split_sequence(seq_records, working_dir)

    # list of dictionaries for each protein containing findername : results
    # entries
    results = [
        {fname : [] for fname,finder in finders.items()}
            for i in range(len(infiles))
    ]

    result_lock = threading.Lock()

    # A queue for all jobs we have to run
    job_queue = queue.Queue()

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

    # launch worker daemon thread
    for i in range(num_threads):
        t = threading.Thread(
            target=finder_worker,
            args=(working_dir, job_queue, results, result_lock)
        )
        t.daemon = True
        t.start()

    logger.debug("Threads active after launching workers: %d",
        threading.active_count())


    logger.info("Processing %d input files in %d jobs.",
        len(infiles), len(joblist))

    # put joblist into job queue, thus starting actual work
    for job in joblist:
        job_queue.put(job)


    job_queue.join()
    logger.info("All jobs returned.")

    #shutil.rmtree(os.path.join(EXECROOT, '..', 'spielwiese'))
    #shutil.copytree(working_dir, os.path.join(EXECROOT, '..', 'spielwiese'))

    return results


######## SET OPEN CONFIGS #########

def set_hhrepid_config_open():
    # construct open configuration for HHrepid
    config = FinderHHrepID.Configuration()
    config.boolopts["-nofilt"] = True
    config.valopts["-P"] = 0.6  # <float> max p-value of suboptimal alignments in all search rounds but the last one (def=0.1)
    config.valopts["-T"] = 0.2  # <float> max total repeat p-value (def=0.001)
    config.valopts["-lmin"] = 0  # [0,inf[  minimal length of repeats to be identified (def=7)
    logger.debug("%s config tokens: %s", finders['hhrepid'].displayname,
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
    logger.debug("%s config tokens: %s", finders['t-reks'].displayname,
                        ", ".join(finders['t-reks'].config.tokens()))
    return config

def set_treks_config_open():
    # construct open configuration for T-Reks
    config = FinderTREKS.Configuration()
    config.valopts["-similarity"] = 0.2
    logger.debug("%s config tokens: %s", finders['t-reks'].displayname,
                        ", ".join(finders['t-reks'].config.tokens()))
    return config

def set_trf_config_open():
    # construct open configuration for TRF
    config = FinderTRF.Configuration()
    config.valopts["Minscore"] = 20
    logger.debug("%s config tokens: %s", finders['trf'].displayname,
                        ", ".join(finders['trf'].config.tokens()))
    return config

def set_trust_config_open():
    # construct open configuration for Trust
    config = FinderTrust.Configuration()
    config.boolopts['-force'] = True
    logger.debug("%s config tokens: %s", finders['trust'].displayname,
                        ", ".join(finders['trust'].config.tokens()))
    return config


def set_xstream_config_open():
    # construct open configuration for XStream
    config = FinderXStream.Configuration()
    #config.boolopts[optname] = True
    config.valopts["-I"] = 0.2
    config.valopts["-i"] = 0.1
    config.valopts["-g"] = 100
    logger.debug("%s config tokens: %s", finders['xstream'].displayname,
                        ", ".join(finders['xstream'].config.tokens()))
    return config


######## RUN A SET OF TRD #########


def run_TRD(sequence_records, sequence_type = 'AA', lFinders = None, default = True):

    ''' Run TRD on sequence_records and return the predicted repeats for each sequence_record and for each tandem repeat detector.
        '''

    # Create temporary working dir
    working_dir = tempfile.mkdtemp()
    logger.debug("repeat_detection_run.run_TRD: Created tempfile: %s", working_dir)

    # Initialise Finders
    Finders(sequence_type, lFinders)

    ## Adjust TRD parameters:
    if not default:
        if sequence_type == 'AA':
            finders['hhrepid'].config = set_hhrepid_config_open()
            finders['trust'].config = set_trust_config_open()
        else:
            finders['trf'].config = set_trf_config_open()
            finders['phobos'].config = set_phobos_config_open()
        finders['t-reks'].config = set_treks_config_open()
        finders['xstream'].config = set_xstream_config_open()

    if sequence_type == 'DNA':
        finders['t-reks'].config = set_treks_config_DNA()

    # start finder. <repeat_detections> has type list(dict("TRD": list(TR)))
    predicted_repeats = run_finders(sequence_records, working_dir=working_dir, num_threads=1)  # [0]  if only the first element is of interest

    #shutil.rmtree(os.path.join(EXECROOT, '..', 'spielwiese'))
    #shutil.copytree(working_dir, os.path.join(EXECROOT, '..', 'spielwiese'))

    # delete temporary directory
    try:
        shutil.rmtree(working_dir)
    except: # I guess the error type is known, and you could be more precise :)
        logging.error("Unexpected error: {0}".format(sys.exc_info()[0]))
        raise

    return predicted_repeats


######## HARDCODED PARAMETERS #########


FINDER_DEFAULT = { "AA": ["HHrepID", "TREKS", "TRUST", "XSTREAM"],
               "DNA": ["Phobos", "TRED", "TREKS", "TRF", "XSTREAM"]
            }

FINDER_FUNCTION_LIST = { "HHrepID": FinderHHrepID(),
                "Phobos": FinderPhobos(),
                "TRED": FinderTRED(),
                "TREKS": FinderTREKS(),
                "TRF": FinderTRF(),
                "TRUST": FinderTrust(),
                "XSTREAM": FinderXStream()
                }

FINDER_LIST = { "HHrepID": FinderHHrepID.name,
                "Phobos": FinderPhobos.name,
                "TRED": FinderTRED.name,
                "TREKS": FinderTREKS.name,
                "TRF": FinderTRF.name,
                "TRUST": FinderTrust.name,
                "XSTREAM": FinderXStream.name
                }