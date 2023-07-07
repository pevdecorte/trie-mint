#!/usr/bin/python3
# Produces the names file from the .ss file.
#
# If a .fa file is provided with the -f switch,
# simply produces a list of fragments,
# but does not include their names.
#
# Sample usages:
#
#   python3 naming.py \
#       data/mm10/mm10-tRNAs-confidence-set.ss \
#       >gen/mm10/mm10-tRNAs.names
#
#   python3 naming.py -f \
#        data/hg19/tRNAspace.Spliced.Sequences.MINTmap_v1.fa \
#       >gen/hg19/hg19-tRNAs.names

import sys
import trnapy
from trnapy import FileFormat


# Generate all tRNA fragments having lengths between range_lower
# and range_upper.
range_lower = 16
range_upper = 50

if __name__ == "__main__":
    usage_message = \
            "\nUsage:\tpython3 postprocess.py [-s] <.ss file>\n" + \
                "\tpython3 postprocess.py -f <.fa file>"
    if len(sys.argv) < 2:
        raise Exception(usage_message)

    file_format = FileFormat.SS
    if sys.argv[1] == "-f" or sys.argv[1] == "-s" or sys.argv[1] == "-t":
        if len(sys.argv) < 3:
            raise Exception(usage_message)
        data_file = open(sys.argv[2], 'r')
        if sys.argv[1] == "-f":
            file_format = FileFormat.FA
        if sys.argv[1] == "-t":
            file_format = FileFormat.TRNA
    else:
        data_file = open(sys.argv[1], 'r')

    if file_format == FileFormat.SS:
        trna_records = trnapy.ReadTRNARecordsFromSSFile(data_file)
    elif file_format == FileFormat.TRNA:
        trna_records = trnapy.ReadTRNARecordsFromTRNAFile(data_file)
    elif file_format == FileFormat.FA:
        trna_records = trnapy.ReadTRNARecordsFromFAFile(data_file)
    data_file.close()

    trna_records.AddCCAToAll()
    trna_records.ExpandAll()

    fragment_sequences_and_types = {}
    def AddFragmentType(fragment):
        if fragment.sequence() not in \
                fragment_sequences_and_types:
            fragment_sequences_and_types[fragment.sequence()] = set()
        if file_format == FileFormat.SS:
            fragment_sequences_and_types[fragment.sequence()].add(
                            fragment.fragment_type())

    for trna in trna_records:
        # If [ACTG] has been prepended.
        if trna.expanded:
            for length in range(range_lower, 1 + min(
                    range_upper, len(trna))):
                fragment = trna.fragment(0, length)
                AddFragmentType(fragment)
        else:
            for start_index in range(len(trna) - range_lower + 1):
                for length in range(range_lower, 1 + min(
                            range_upper, len(trna) - start_index)):
                    fragment = trna.fragment(start_index, length)
                    AddFragmentType(fragment)

    for sequence, types in fragment_sequences_and_types.items():
        if file_format == FileFormat.SS:
            print(sequence + "\t", end="")
            print(*types, sep=", ")
        else:
            print(sequence)

