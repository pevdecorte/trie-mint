#!/usr/bin/python3
# Usage:
#   python3 reduced_postprocess.py <.ss file>
#       <.names file> [MATCHES FILES]
#
# The lines in the .matches files should have the following
# format.
#
#       <fragment> <position in chromosome>
#           <position in trna> [!]<original indices in trna>
#
#   In addition, the file should have a header indicating
#   the name of the haystack file.
#
# For each <fragment> <name> line in the .names file, we 
# print the following line.
#
# <fragment> <Y|N>
#
# If the fragment is exclusive to tRNA space, the fragment
# is followed by a 'Y'. Otherwise it is followed by an 'N'.

import os.path
import re
import sys
import trnapy

from trnapy import FileFormat
from trnapy import InverseComplement


# Returns True iff interval is contained in the union
# of the intervals listed in interval_list.
def ContainedInIntervalList(interval_list, interval):
    l = 0; r = len(interval_list) - 1
    while r - l > 1:
        m = (l + r) // 2
        if interval_list[m][0] < interval[0]:
            l = m
        elif interval_list[m][0] >= interval[0]:
            r = m

    if interval_list[l][0] <= interval[0] < interval[1] <= \
            interval_list[l][1]:
        return True
    if interval_list[r][0] <= interval[0] < interval[1] <= \
            interval_list[r][1]:
        return True
    if interval_list[l][0] <= interval[0] <= interval_list[l][1] \
            and interval_list[r][0] - 1 == interval_list[l][1] \
            and interval_list[r][0] <= interval[1] <= interval_list[r][1]:
        return True

    return False


# Returns True iff the interval in match_dict_key is contained
# in tRNA space, as described by trna_space.
def ContainedInTRNASpace(trna_space, match_dict_key):
    trna_space_key = (match_dict_key[0], match_dict_key[1])
    interval = match_dict_key[2]
    if trna_space_key not in trna_space:
        return False
    space = trna_space[trna_space_key]
    
    return ContainedInIntervalList(space, interval)

# Parses a tRNA-unaware match line of the form:
# GGGGTTGGGGATTTAG 31653 4 0-81
# Returns the sign and the subinterval of original nucleotides
# (i.e. CCA and possibly the leading extra nucleotide are not counted.
# The returned value has the form
# (fragement_sequence, sign, interval_start, interval_end), 
def ParseTRNAUnawareMatch(match_line):
    seq, pos_in_genome, pos_in_trna, original_interval = \
            match_line.strip().split(' ')
    pos_in_genome = int(pos_in_genome)
    pos_in_trna = int(pos_in_trna)

    if original_interval[0] == '!':
        sign = '-'
        original_interval = original_interval[1:]
    else:
        sign = '+'
    original_interval_start, original_interval_end = \
            original_interval.strip().split('-')
    original_interval_start = int(original_interval_start)
    original_interval_end = int(original_interval_end)

    interval_start = pos_in_genome
    start_diff = original_interval_start - pos_in_trna
    if start_diff > 0:
        interval_start += start_diff

    interval_end = pos_in_genome + len(seq) - 1
    end_diff = pos_in_trna + len(seq) - 1 - original_interval_end
    if end_diff > 0:
        interval_end -= end_diff

    return (seq, sign, interval_start, interval_end)


if __name__ == "__main__":
    usage_message = "Usage:\tpython3 reduced_postprocess.py [-fs] " +\
            "[<.ss file>|<.fa file>|<.trna file>] <.names file> [.matches FILES]"

    if len(sys.argv) < 3:
        raise Exception(usage_message)

    file_format = FileFormat.SS
    switch_provided = False
    if sys.argv[1] == "-f" or sys.argv[1] == "-s" or sys.argv[1] == "-t":
        switch_provided = True
        if len(sys.argv) < 3:
            raise Exception(usage_message)
        records_file = open(sys.argv[2], 'r')
        if sys.argv[1] == "-f":
            file_format = FileFormat.FA
        elif sys.argv[1] == "-t":
            file_format = FileFormat.TRNA
    else:
        records_file = open(sys.argv[1], 'r')

    if file_format == FileFormat.SS:
        trna_records = trnapy.ReadTRNARecordsFromSSFile(records_file)
    elif file_format == FileFormat.FA:
        trna_records = trnapy.ReadTRNARecordsFromFAFile(records_file)
    elif file_format == FileFormat.TRNA:
        trna_records = trnapy.ReadTRNARecordsFromTRNAFile(records_file)

    trna_records.AddCCAToAll()
    trna_records.ExpandAll()

    # Matches against fragments of virtual tRNA's correspond
    # to matches on the negative strand. Virtual tRNA's are
    # not used in defining tRNA space.
    virtual_trna_records = trna_records.InverseComplements()

    trna_space = trnapy.ConstructTRNASpace(trna_records)

    outside_trna_space = set()
    for matches_filename in sys.argv[4:] if switch_provided else \
            sys.argv[3:]:

        with open(matches_filename, 'r') as matches_file:
            # Parses out the chromosome name from the haystack
            # file name.
            chromosome_filename = next(matches_file).strip()
            chromosome_filename = os.path.basename(chromosome_filename)
            m = re.search("^[^.]*", chromosome_filename)
            chromosome = m.group(0)

            for line in matches_file:
                seq, sign, interval_start, interval_end = \
                        ParseTRNAUnawareMatch(line)

                match_dict_key = (chromosome, sign,
                    (interval_start, interval_end))

                if not ContainedInTRNASpace(trna_space, match_dict_key):
                    outside_trna_space.add(seq if sign == "+" else InverseComplement(seq))

    with open(sys.argv[3] if switch_provided else sys.argv[2],
            'r') as names_file:
        for line in names_file:
            fields = line.strip().split()
            if fields[0] in outside_trna_space:
                print(fields[0], "N", sep='\t')
            else:
                print(fields[0], "Y", sep='\t')

