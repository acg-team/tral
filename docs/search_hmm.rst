.. _search_hmm:

Searching for particular repeats
********************************

Overview
========

The search package provides python modules and command line tools for identifying
tandem repeats in sequence databases. Typically processing consists of first
running ``search_hmm`` to identify instances of a query repeat
in a sequence database. This is relatively slow (hours to days
against Uniprot) but can be trivially parallelized to run on a cluster.

After scoring the hmm against every sequence in the database, ``filter_hmm`` is
used to identify significant results.

search_hmm
==========

The search module can be run directly:
::

    python -m tral.search.search_hmm -h

The query HMM should be in HMMR format. The ``hmmbuild`` tool may be helpful for
generating HMM files from multiple alignments. The HMM should represent a single
repeat. Circularly permuted edges are automatically inserted to permit multiple
sequential repeats.

The database should be in FASTA format. It may be gzip compressed.

Results are produced in TSV format, and include

* Sequence identifier (header line from FASTA file)
* Viterbi probability (probability that this sequence was emitted by the query
  HMM)
* Log-odds ratio (viterbi probability divided by the expected probability for
  a sequence of this length
* HMM states (Maximum likelihood sequence of states through the HMM, as
  determined by the Viiterbi algorithm

Example:
::

    python -m tral.search.search_hmm query.hmm database.faa.gz results.tsv

When running on a cluster, it may be useful to search just a subset of the
database. For instance, here's a simplified submission script using the ``-n``
and ``-s`` parameters in an array job:
::

    #!/bin/bash
    #BSUB -J tral[1-10]
    n=1000
    start=$(( (${LSB_JOBINDEX:-1}-1) * $n ))
    python -m tral.search.tral_search \
        -n $n -s $start \
        query.hmm database.fasta.gz results.${LSB_JOBINDEX:-1}.tsv

Choosing an appropriate log-odds ratio depends on the database size, amino
acid content, HMM length, and HMM composition. To help calibrate this, the ``-r``
option is available to automatically shuffle the input sequences. The results
then represent all negative hits, allowing the False Discovery Rate to be
estimated for various thresholds in the filtering step.

Full usage
----------
::

    usage: search_hmm.py [-h] [-v] [-s START] [-n DBSIZE] [-r]
                        hmm database results

    Align a cpHMM against a database using TRAL

    positional arguments:
    hmm                   cpHMM filename
    database              database to search, in fasta format (may be gzip
                            compressed)
    results               results file

    optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         Long messages
    -s START, --start START
                            Index of database sequence to start with (default: 0)
    -n DBSIZE, --dbsize DBSIZE
                            Number of database sequences to include (default: 0 to
                            include all)
    -r, --shuffle         Shuffle database sequences



filter_hmm
==========

The results from ``search_hmm`` include both significant and insignificant hits.
Use ``filter_hmm`` to filter the results based on length and log-odds criteria.
::

    python -m tral.search.filter_hmm -h

**Output**. The tool can produce three kinds of output:

* ``-f`` **FASTA**: Subset of the input database corresponding to significant hits
* ``-o`` **TSV**: Subset of the ``search_hmm`` results corresponding to
  significant hits
* ``-x`` **TREKS**: All significant repeats in T-REKS format.

The TREKS format is most useful for exporting results to other programs.

**Identifiers**. By default, the identifier is extracted from the database
header lines (e.g. ``>tr|id|description...``). For compatibility with other
programs, the ``--preserve-header`` option can also be used to match sequences
based on the full header.

**Thresholds**. Thresholds can be set for either:

* ``-r`` **Number of repeats** (float). Repeats are computed based on the HMM match
  states. This means that partial repeats are counted as a fraction of the
  query length.
* ``-t`` **Log-odds threshold** (float). Results with log-odds less than this value
  are filtered out.

**Memory usage**. Memory usage should be low if a single output is used. If
several output options are selected then the significant hits are stored in
memory. If this is problematic, ``filter_hmm`` can be run separately for
each desired output.


**Example:**
::

    python -m tral.search.filter_hmm \
        -f filtered.faa \
        -o filtered.tsv \
        -x filtered.tdr \
        -r 2 -t 8.0 \
        results.tsv database.faa.gz


Full usage
----------
::

    usage: filter_hmm.py [-h] [-f FILTERED_FASTA] [-o FILTERED_TSV]
                        [-x FILTERED_TREKS] [-r MIN_REPEATS] [-t LOG_ODDS]
                        [--preserve-header] [-v]
                        hits database

    Filter hits from search_hmm

    positional arguments:
    hits                  TSV file, as produced by search_hmm, containing HMM
                            hits
    database              database to filter, in fasta format (may be gzip
                            compressed)

    optional arguments:
    -h, --help            show this help message and exit
    --preserve-header     Include the full header in FASTA output. Otherwise,
                            just the identifier is used to match TSV and TREKS.
    -v, --verbose         Long messages

    outputs:
    -f FILTERED_FASTA, --filtered-fasta FILTERED_FASTA
                            output fasta file, filtered by hits
    -o FILTERED_TSV, --filtered-tsv FILTERED_TSV
                            Filtered TSV file
    -x FILTERED_TREKS, --filtered-treks FILTERED_TREKS
                            Filtered T-Reks file

    Thresholds:
    -r MIN_REPEATS, --min-repeats MIN_REPEATS
                            Minimum number of repeats
    -t LOG_ODDS, --log-odds LOG_ODDS
                            Threshold for minimum log-odds ratio


Python API
==========

See also the :ref:`Search code documentation <tral/search>`
