# (C) 2011, Alexander Korsunsky
# (C) 2011-2015 Elke Schaper

"""
    :synopsis: Input/output for repeat detection algorithms

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>

"""

import logging
import re
import csv

LOG = logging.getLogger(__name__)


class RepeatRegion:

    def __init__(self, protein_id="", begin=None, msa=None):
        self.protein_id = protein_id
        self.begin = begin
        if msa is None:
            msa = []
        self.msa = msa


def tred_get_repeats(infile):
    """ Read repeats from a TRED standard output (stdout) file stream successively.

    Read repeats from a TRED standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of TRED output file::

         Start: start End: \d+ Length: \d+

        ( \d repeat_unit \d
        ( alignment_indicator )?)*

    Args:
        infile (file stream): File stream of output1 from
            tred1 myDNA.faa intermediate_output
            tred2 myDNA.faa intermediate_output output1 output2

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRED output syntax.
    """

    # pattern for start index
    pat_start = re.compile(r"Start: (\d+)")
    pat_repeat_unit = re.compile(r"\d+\s+([ACGT-]+)")
    pat_alignment_indicator = re.compile(r"\s+([ACGT-]+)")

    # Our possible parser states:
    #
    # 1: searching for start
    # 2: searching for first repeat unit
    # 3: searching for all further repeat units

    identifier = ''
    state = 1
    repeat_units = []
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])
        if 1 == state:
            region = RepeatRegion(protein_id=identifier, msa=[])
            match = pat_start.match(line)
            if match:
                LOG.debug(" * (1->2) Found start")
                LOG.debug("Start: %s", match.group(1))
                region.begin = int(match.group(1))
                state = 2

        if 2 == state:
            match = pat_repeat_unit.match(line)
            if match:
                LOG.debug(" * (2->3) Found first repeat_unit")
                repeat_units.append((match.group(1), True))
                state = 3
                continue

        if 3 == state:
            match = pat_repeat_unit.match(line)
            if match:
                LOG.debug(" * (3->3) Found another repeat_unit")
                repeat_units.append((match.group(1), True))
            else:
                match = pat_alignment_indicator.match(line)
                if match:
                    LOG.debug(" * (3->3) Found an alignment_indicator unit")
                    repeat_units.append((match.group(1), False))
                else:
                    LOG.debug(" * (3->1) Found end of repeat (yielding)")
                    state = 1
                    region.msa = tred_msa_from_pairwise(repeat_units)
                    yield region


def tred_msa_from_pairwise(repeat_units):
    """ Construct a MSA from pairwise alignments.

    Construct a MSA from pairwise alignments. At the moment, gaps following the repeat are
    not added. However, these gaps are added automatically when a ``Repeat`` instance is
    created.

    Args:
        repeat_units (list of str): Read in from TRED output files

    Returns:
         (list of str)

    .. todo:: Is the Args format mentioned correctly?
    """

    pat_gap = re.compile(r"(-*)")
    result = []
    index = 0
    for iR in range(len(repeat_units)):
        ru = repeat_units[iR]

        # The next repeat unit
        if ru[1] == True:
            result.append('-' * index + ru[0])
            # How many gaps in the beginning of this repeat unit?
            index += len(pat_gap.match(ru[0]).group())

        # The next alignment indicator
        else:
            for iL in range(len(ru[0])):
                if ru[0][iL] == '-':
                    # enter a gap between
                    # the index + iL and index + iL + 1 character
                    # in each repeat unit in result so far:
                    result = [iRU[:index + iL] + '-' + iRU[index + iL:]
                              for iRU in result]

    return result


def treks_get_repeats(infile):
    """ Read repeats from a T-REKS standard output (stdout) file stream successively.

    Read repeats from a T-REKS standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of T-REKS output file::

        protein ::=
            ">" identifier
            repeat*
        #
        repeat ::=
            repeat_header
            sequence*
            "*"+
        #
        repeat_header ::= "Length:" integer "residues - nb:" integer  "from"  integer "to" integer "- Psim:"float "region Length:"integer


    Args:
        infile (file stream): File stream to the file of the standard output of T-Reks

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout T-REKS output syntax.
    """

    # pattern for protein identifier
    pat_identifier = re.compile(r">(\S+)")

    # pattern for repeat properties
    pat_repeat_header = re.compile(
        r"Length: \d+ residues - nb: (\d+)  from  (\d+) to (\d+) - Psim:([\d\.]+) region Length:(\d+)")

    pat_repeat_end = re.compile(r"\*+")

    # pattern for repeat sequence
    # FIXME stuff that occurs here is not just A-Z but a set of possible amino
    # acid symbols
    pat_sequence = re.compile(r"([A-Z-]+)")

    # Our possible parser states:
    #
    # 1: state between repeats
    #   entry: reset repeat
    #   expect repeat_header(goto 2) OR identifier(store identifier, goto 1)
    # 2: state for multiple sequence alignment line
    # expect sequence(goto 3, append sequence) OR repeat_end(return repeat,
    # goto 1)

    state = 1
    identifier = ""

    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        if 1 == state:
            region = RepeatRegion(protein_id=identifier, msa=[])
            LOG.debug("msa: %s", "\n".join(region.msa))
            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" * (1->2) Found properties")
                region.begin = int(match.group(2))
                state = 2

            match = pat_identifier.match(line)
            if match:
                LOG.debug(" * (1->1) Found identifier")
                identifier = match.group(1)
                state = 1

        elif 2 == state:
            match = pat_sequence.match(line)
            if match:
                LOG.debug(" * (2->2) Found MSA line (appending)")
                region.msa.append(match.group(1))
                state = 2
            match = pat_repeat_end.match(line)
            if match:
                LOG.debug(" * (2->1) Found end of repeat (yielding)")
                state = 1
                yield region

        else:
            raise AssertionError("Huh? Unknown parser state " +
                                 str(state))


def xstream_get_repeats(infile):
    """ Read repeats from a XSTREAM output xls chart

    Read repeats from a XSTREAM output xls chart

    Postcondition: infile points to EOF.

    Args:
        infile (file stream): File stream to read XSTREAM output from a xls chart.

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.
    """

    # The infile is luckily enough in csv format:
    reader = csv.reader(infile, dialect='excel-tab')

    # read first row with the fieldnames
    header = next(reader)

    LOG.debug("Header has %d fields.", len(header))

    if len(header) != 8 and len(header) != 10:
        raise ValueError(
            "XStream output file seems to be malformed. "
            "Make sure to open this function with the generated xls chart.")

    # when XSTREAM is fed files with multiple sequences, first field is "identifier",
    # otherwise it is "seq length"
    if header[0] == "identifier":
        field_offset = 2
    else:
        field_offset = 0

    for row in reader:
        region = RepeatRegion()
        if field_offset:
            region.protein_id = row[0]

        region.begin = int(row[0 + field_offset])
        region.msa = row[4 + field_offset].split()
        LOG.debug("Found repeat, yielding.")
        if len(region.msa) >= 2:
            yield region


def trust_fill_repeats(msa, begin, sequence, maximal_gap_length=20):
    ''' return a trust msa that has no longer indels than maximal_gap_length,
    that contains the indel characters even when not part of the trust output file.
    Background trust returns tandem repeats, but also distant repeats. '''
    gapless_msa = [repeat_unit.replace('-', '').upper() for repeat_unit in msa]
    sequence = sequence.upper()

    # Find the start and end positions of the predicted repeat units
    position = [(begin - 1, begin + len(gapless_msa[0]) - 2)]
    for repeat_unit in gapless_msa[1:]:
        find_index = sequence[position[-1][1] + 1:].find(repeat_unit)
        # Repeat unit could not be found in sequence -> Discard the repeat
        if find_index == -1:
            return None, None
        repeat_unit_begin = find_index + position[-1][1] + 1
        position.append(
            (repeat_unit_begin, repeat_unit_begin + len(repeat_unit) - 1))

    # Derive the start and end positions of the gaps
    gap_position = [(i[1] + 1, j[0] - 1)
                    for i, j in zip(position[:-1], position[1:])]
    gaps = [i[1] - i[0] + 1 for i in gap_position]

    # Filter out repeat units that are further apart than maximal_gap_length
    gap_valid = ''.join(
        ['1' if i_gap <= maximal_gap_length else '0' for i_gap in gaps])

    count_valid_pairs = [len(m.group())
                         for m in re.finditer(re.compile('1+'), gap_valid)]
    # All repeat units are further apart than maximal_gap_length? -> Discard
    # the repeat
    if len(count_valid_pairs) == 0:
        return None, None

    # Choose the sequence of pairs closer then maximal_gap_length that is
    # longest
    valid_index = gap_valid.find('1' * max(count_valid_pairs))

    # Shorten the predicted msa accordingly.
    msa = msa[valid_index:valid_index + max(count_valid_pairs) + 1]
    gaps = gaps[valid_index:valid_index + max(count_valid_pairs) + 1]
    gap_position = gap_position[
        valid_index:valid_index +
        max(count_valid_pairs) +
        1]
    position = position[valid_index:valid_index + max(count_valid_pairs) + 1]

    # Add missing sequence to the repeat units
    repeat_unit_length = len(msa[0])
    gap_count_before = 0
    for i, i_gap in enumerate(gaps):
        gap_count_after = gap_count_before + i_gap
        msa[i] += '-' * (gap_count_before) + sequence[gap_position[i][0]:gap_position[i][1] + 1] + '-' * (sum(gaps) - gap_count_after)
        gap_count_before = gap_count_after
    msa[-1] += '-' * sum(gaps)

    return msa, position[0][0] + 1


def trust_get_repeats(infile):
    """ Read repeats from a TRUST standard output (stdout) file stream successively.

    Read repeats from a TRUST standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of TRUST standard output::

        protein ::=
            ">" identifier
            (repeat_types)*
            "//"
        #
        repeat_types ::=
            "REPEAT_TYPE" integer
            "REPEAT_LENGTH" integer
            (repeat_info)*
            (">Repeat " integer
            sequence)*
        #
        repeat_info ::=
            integer integer [integer] [integer]

    Args:
        infile (file stream): File stream from TRUST standard output.

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRUST output syntax.
    """

    pat_identifier = re.compile(r">(\S+)")
    pat_repeat_type = re.compile(r"REPEAT_TYPE \d+")
    pat_repeat_length = re.compile(r"REPEAT_LENGTH \d+")
    pat_repeat_info = re.compile(r"(\d+) (\d+).*$")
    pat_repeat_header = re.compile(r">Repeat \d+")
    # FIXME find proper character set here
    pat_repeat_sequence = re.compile("([A-Za-z-]+)")
    pat_protein_end = re.compile("//")

    # Our possible parser states:
    #
    # 0: initial state
    #   expect: identifier(store identifier, goto 1)
    # 1: before the beginning of a repeat
    #   expect: "REPEAT_TYPE"(goto 2) or "//"(goto 0)
    # 2: state after having found "REPEAT_TYPE"
    #   entry: reset state of return value (repeat)
    #   expect: "REPEAT_LENGTH"(goto 3)
    # 3: first line of repeat info
    #   expect: repeat_info(store begin, goto 4) or ">Repeat"(goto 5)
    # 4: continuing to read repeat info
    #   expect: repeat_info(goto 4) or ">Repeat"(goto 5)
    # 5: part of MSA sequence
    #   expect: sequence(goto 6)
    # 6: After first sequence
    # expect: ">Repeat"(goto 5) or "REPEAT_TYPE"(return repeat, goto 2) or
    # "//"(return repeat, goto 0)

    state = 0
    region = RepeatRegion()
    identifier = ""

    def strip_comments(line):
        pat_comment = re.compile(r"\s*#.*$")
        return pat_comment.sub("", line)

    for i, line in enumerate(infile):
        line = strip_comments(line)
        if not line or line == '\n':
            continue

        LOG.debug("Line %d: %s", i, line[0:-1])

        if 0 == state:
            match = pat_identifier.match(line)
            if match:
                LOG.debug(
                    " *(0->1) Found identifier (storing \"%s\")",
                    match.group(1))

                identifier = match.group(1)
                state = 1
                continue

        elif 1 == state:
            match = pat_repeat_type.match(line)
            if match:
                LOG.debug(" *(1->2) Found REPEAT_TYPE")
                state = 2
                continue

            match = pat_protein_end.match(line)
            if match:
                LOG.debug(" *(1->0) Found protein end")
                state = 0
                continue

        elif 2 == state:
            # Entry Action: clear msa cache
            region = RepeatRegion(protein_id=identifier)

            match = pat_repeat_length.match(line)
            if match:
                LOG.debug(" *(2->3) Found REPEAT_LENGTH")
                state = 3
                continue

        elif 3 == state:
            match = pat_repeat_info.match(line)
            if match:
                LOG.debug(
                    " *(3->4) Found repeat_info (storing begin: \"%s\")",
                    match.group(1)
                )

                region.begin = int(match.group(1))
                state = 4
                continue

            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(3->5) Found repeat_header")
                state = 5
                continue

        elif 4 == state:
            match = pat_repeat_info.match(line)
            if match:
                LOG.debug(" *(4->4) Found repeat_info")
                state = 4
                continue

            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(4->5) Found repeat_header")
                state = 5
                continue

        elif 5 == state:
            match = pat_repeat_sequence.match(line)
            if match:
                LOG.debug(
                    " *(5->6) Found sequence (storing \"%s\")",
                    match.group(1))

                region.msa.append(match.group(1))

                state = 6
                continue

        elif 6 == state:
            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(6->5) Found repeat_header")
                state = 5
                continue
            match = pat_repeat_type.match(line)
            if match:
                LOG.debug(" *(6->2) Found REPEAT_TYPE, yielding.")

                state = 2
                if len(region.msa) >= 2:
                    yield region
                continue
            match = pat_protein_end.match(line)
            if match:
                LOG.debug(" *(6->0) Found protein end, yielding.")
                state = 0
                if len(region.msa) >= 2:
                    yield region
                continue

        else:
            raise AssertionError("Huh? Unknown parser state " + str(state))

################################## TRF - Benson ##########################


def trf_get_repeats(infile):
    """ Read repeats from a TRF txt.html file stream file stream successively.

    Read repeats from a TRF txt.html file stream file stream successively.
    Postcondition: infile points to EOF.

    TRF output file syntax::

        Sequence: ``identifier``
             Indices: ``begin``--``end``
             \d [a-zA-Z]+
        #
             begin (repeat)*
             1  (consensus)*
        #
          (( \d (repeat)*
             \d  (consensus)*
          )?
             \d (repeat)*
             1  (consensus)*
          )+
             \d [a-zA-Z]+
        #
            ``Statistics``

    Args:
        infile (file stream): File stream from TRF output txt.html.
    (generated via e.g. ./trf404.linux64.exe FASTAFILE 2 7 7 80 10 50 500 -d > /dev/null
    If the -h flag is set, no .txt.html output is produced)

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRF output syntax.
    .. todo:: Does not search for the sequence identifier at current!
    """

    # find the name of the sequence ## CURRENTLY NOT IMPLEMENTED
    #pat_identifier = re.compile("Sequence: (\S+)")

    # find coordinates of the repeat region in the protein
    pat_coordinates = re.compile("Indices: (\d+)--(\d+)")

    # find a part of a repeat unit and its coordinate
    pattern_seq = re.compile("\s+(\d+) ([ACGT\- ]+)")

    # find the final tag 'Statistics'
    pat_statistics = re.compile("Statistics")

    # Our possible parser states:
    #
    # state 0: searching for identifier -> 1  # not necessary when sequence identifier is known
    # state 1: searching for repeat region coordinates -> 2
    # state 2: searching for beginning of MSA & save sequence to tmpMSA-> 4
    # state 3: new repeat unit: save sequence to tmpMSA -> 6
    # state 4: new repeat unit: save sequence to tmp_consensus -> 5
    # state 5:
    #           if sequence: save sequence to tmpMSA -> 6
    #           if end: save tmpMSA to preMSA; save tmp_consensus to consensus; Yield repeat region -> 1
    # state 6: check: new repeat unit? save tmp_consensus
    #  1: use last tmpMSA entry for new tmpMSA;
    #     save tmpMSA to preMSA;
    #     save tmp_consensus to consensus;
    #     save sequence to new tmp_consensus -> 3
    #  0: save sequence to tmp_consensus -> 5

    state = 1
    #identifier = ""  # Currently not implemented.
    preMSA = []
    consensus = []
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        # CURRENTLY NOT IMPLEMENTED
        # if state == 0: #searching for sequenceMSA identifier
        #  tmp = pat_identifier.search(line)
        #  if tmp:
        #    state = 1
        #    identifier = tmp.group()
        #  continue
        # elif state == 1: #searching for TR boundaries (indices)

        if 1 == state:  # searching for repeat region coordinates
            search = pat_coordinates.search(line)
            if search:
                LOG.debug(" * (1->2) Found coordinates")
                state = 2
                region = RepeatRegion()
                region.begin = int(search.group(1))
                #region_end = search.group(2)
                short = False

        # searching for beginning of MSA & save sequence to tmpMSA-> 4
        elif state == 2:
            match = pattern_seq.match(line)
            if match and match.group(1) == str(region.begin):
                LOG.debug(" *(2->4) Found first row of first MSA repeat unit")
                state = 4
                if len(match.group(2).strip().split(" ")) > 1:
                    short = True
                    preMSA = match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA = [match.group(2).strip().split(" ")[0]]
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))

        elif state == 3:  # new repeat unit: save sequence to tmpMSA -> 4
            match = pattern_seq.match(line)
            if match:
                LOG.debug(" *(3->5) Found first row of new repeat unit")
                state = 6
                if short:
                    preMSA += match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA.append(match.group(2).strip().split(" ")[0])
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))
            # if end: save tmpMSA to preMSA; save tmp_consensus to consensus;
            # Yield repeat region -> 1
            if pat_statistics.search(line):
                LOG.debug(
                    " *(5->1) Encountered 'Statistics': No more repeats, yielding.")
                state = 1
                if not short:
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []

        elif state == 4:  # new repeat unit: save sequence to tmp_consensus -> 5
            match = pattern_seq.match(line)
            if match:
                LOG.debug(
                    " *(4->5) Found first consensus row of the repeat unit")
                state = 5
                if short:
                    consensus = match.group(2).strip().split(" ")
                    LOG.debug(
                        "Repeat unit is short;  consensus: %s",
                        str(consensus))
                else:
                    tmp_consensus = [match.group(2).strip().split(" ")[0]]
                    LOG.debug(" tmp_consensus: %s", str(tmp_consensus))

        elif state == 5:  # SEARCHING FOR MSA ROW
            # if sequence: save sequence to tmpMSA -> 6
            match = pattern_seq.match(line)
            if match:
                LOG.debug(" *(5->6) Found a MSA row")
                state = 6
                if short:
                    preMSA += match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA.append(match.group(2).strip().split(" ")[0])
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))

            # if end: save tmpMSA to preMSA; save tmp_consensus to consensus;
            # Yield repeat region -> 1
            if pat_statistics.search(line):
                LOG.debug(
                    " *(5->1) Encountered 'Statistics': No more repeats, yielding.")
                state = 1
                if not short:
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []

        elif state == 6:  # new repeat unit? ## SEARCHING FOR CONSENSUS ROW
            match = pattern_seq.match(line)
            # 1: save tmp_consensus -> 3
            if match and match.group(1) == '1':
                LOG.debug(
                    " *(6->3) Found a consensus row of a new repeat unit")
                state = 3
                if short:  # NEEDS TO BE CODED
                    consensus += match.group(2).strip().split(" ")
                    LOG.debug(
                        "Repeat unit is short; consensus: %s",
                        str(consensus))
                else:
                    newMSA = tmpMSA.pop()
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                    tmpMSA = [newMSA]
                    tmp_consensus = [match.group(2).strip().split(" ")[0]]

            # 0: save sequence to tmp_consensus -> 5
            elif match:
                LOG.debug(" *(6->5) Found a consensus row")
                state = 5
                tmp_consensus.append(match.group(2).strip().split(" ")[0])

            # YIELD
            # aha! there should have been a consensus sequence, but there is
            # not. Hence we are finished with this repeat!
            else:
                LOG.debug(' *(6->1) No consensus row: repeat finished')
                state = 1
                if short:
                    preMSA = preMSA[:-1]
                else:
                    preMSA.append(
                        "".join(
                            tmpMSA[
                                0:-
                                1]))  # The last tmpMSA entry was not a repeat unit
                    consensus.append("".join(tmp_consensus))
                LOG.debug(" preMSA: %s", str(preMSA))
                LOG.debug(" consensus: %s", str(consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []


def getMSA(sequenceMSA, consensusMSA):
    """ Derive the MSA from a strange combination of consensusMSA and sequenceMSA in TRF
    (Benson) txt.html output files

    Args:
        sequenceMSA (?):
        consensusMSA (?):

    Returns:
         msa (list of str): The multiple sequence alignment predicted by TRF.
    """

    msa = [""] * len(sequenceMSA)

    while consensusMSA:
        # CHECK for insertions
        insertion = 1
        while insertion and consensusMSA:
            insertion = 0
            for i_con in consensusMSA:
                if i_con and i_con[0] == "-":
                    insertion = 1
                    break
            # INCLUDE insertions into the msa
            if insertion:
                for i in range(len(consensusMSA)):
                    if consensusMSA[i] and consensusMSA[i][0] == "-":
                        msa[i] += sequenceMSA[i][0]
                        sequenceMSA[i] = sequenceMSA[i][1:]
                        consensusMSA[i] = consensusMSA[i][1:]
                    else:
                        msa[i] += "-"

        # CHECK for deletions and normal sequence
        if not consensusMSA[0]:
            break

        for i in range(len(consensusMSA)):
            # The last repeat unit can be shorter than the ones before
            if not sequenceMSA[i]:
                break

            msa[i] += sequenceMSA[i][0]
            sequenceMSA[i] = sequenceMSA[i][1:]
            consensusMSA[i] = consensusMSA[i][1:]

    return msa

# ############## HHrepID - Soeding ###########################################


def hhpredid_get_repeats(infile):
    """ Read repeats from a HHREPID standard output (stdout) file stream successively.

    Read repeats from a HHREPID standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of HHREPID standard output::

        protein ::=
             begin"-"\d    "+"\d repeatUnit
           ( \d"-"\d    "+"\d repeatUnit )+

    Args:
        infile (file stream): File stream from HHREPID standard output.
        [Generated by e.g.: ./hhrepid_32 -i FASTAFILE -v 0 -d cal.hhm -o INFILE]

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout HHREPID output syntax.
    """

    # find a part of a repeat unit and its first coordinate
    # minus or \minus?

    pattern_repeat_unit_count = re.compile("Repeats\s+(\d+)")
    pattern_seq = re.compile("[A-Z]+(\d+).*(\d+)\-.*\+[\d]+ ([\-a-zA-Z.]+)")

    # Our possible parser states:

    # state1: Find number of repeat units n
    # state2: Find first (partial) row of the MSA
    # state3: Find all other (partial) rows of the MSA

    region = None
    state = 1
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        if 1 == state:  # Find 'Repeats' marker of new repeat
            search = pattern_repeat_unit_count.search(line)
            if search:
                LOG.debug(" *(1->2) Found repeat")
                state = 2
                n = int(search.group(1))

        elif 2 == state:  # Find first (partial) row of the MSA
            search = pattern_seq.search(line)
            if search:
                LOG.debug(" *(2->3) Found first repeat unit (part)")
                state = 3
                region = RepeatRegion()
                region.begin = int(search.group(2))
                region.msa = [""] * n
                region.msa[int(search.group(1)) -
                           1] = search.group(3).replace('.', '-').upper()

        elif 3 == state:  # Find all other (partial) rows of the MSA
            search = pattern_seq.search(line)
            if search:
                LOG.debug(" *(3->3) Found other repeat unit (part)")
                region.msa[int(search.group(1)) -
                           1] += search.group(3).replace('.', '-').upper()
            else:
                search = pattern_repeat_unit_count.search(line)
                if search:
                    LOG.debug(" *(3->2) Yield Repeat, begin next")
                    state = 2
                    n = int(search.group(1))
                    if len(region.msa) >= 2:
                        yield region
                        region = None
                    else:
                        log.warning(
                            "HHPREDID: Msa too short %s", str(
                                region.msa))

    # Yield final repeat region.
    if not region is None:
        if len(region.msa) >= 2:
            yield region
        else:
            log.warning("HHPREDID: Msa too short %s", str(region.msa))

####################################### Phobos TRF  ######################


def phobos_get_repeats(infile):
    """ Read repeats from a PHOBOS output file stream successively.

    Read repeats from a PHOBOS output file stream successively.
    Postcondition: infile points to EOF.

    Args:
        infile (file stream): File stream from PHOBOS output.

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Show PHOBOS output syntax.
    """

    pattern_begin = re.compile("(\d+) :\s+\d")
    pattern_seq = re.compile("([\-ACGT]+)")

    # Our possible parser states:
    #
    # state 1: Find TR begin
    # state 2: Find first repeat unit
    # state 3: Find repeat units

    state = 1
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])
        if 1 == state:  # Find TR offset
            search = pattern_begin.search(line)
            if search and search.groups()[0] != None:
                LOG.debug(" *(1->2) Found tandem repeat begin")
                state = 2
                region = RepeatRegion()
                region.begin = int(search.groups()[0])
                region.msa = []

        elif 2 == state:  # Find all other repeat units
            match = pattern_seq.search(line)
            if match and match.groups()[0] != None:
                LOG.debug(" *(2->3) Found first repeat unit")
                region.msa.append(match.groups()[0])
                state = 3

        elif 3 == state:  # Find all other repeat units
            match = pattern_seq.search(line)
            if match and match.groups()[0] != None:
                LOG.debug(" *(3->3) Found a repeat unit")
                region.msa.append(match.groups()[0])
            else:
                LOG.debug(" *(3->1) repeat region finished, yielding.")
                state = 1
                if len(region.msa) >= 2:
                    yield region
                else:
                    log.warning("phobos: Msa too short %s", str(region.msa))
