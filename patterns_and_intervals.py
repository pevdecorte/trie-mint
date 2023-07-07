#!/usr/bin/python3
# Create a .patterns file from an .ss file or a .fa file.
#
# Use the -f switch when passing a .fa file.
#
# Sample usages:
#   python3 patterns_and_intervals.py -s \
#        data/mm10/mm10-tRNAs-confidence-set.ss
#   python3 patterns_and_intervals.py -f \
#        data/hg19/tRNAspace.Spliced.Sequences.MINTmap_v1.fa
#
# By default the script assumes a .ss file is being passed.

import sys
import trnapy
from trnapy import FileFormat

if __name__ == "__main__":
    usage_message = ".ss file or .fa must be provided"
    if len(sys.argv) < 2:
        raise Exception(usage_message)

    file_format = FileFormat.SS
    if sys.argv[1] == "-f" or sys.argv[1] == "-s" or sys.argv[1] == "-t":
        if len(sys.argv) < 3:
            raise Exception(usage_message)
        data_file = open(sys.argv[2], 'r')
        if sys.argv[1] == "-f":
            file_format = FileFormat.FA
        elif sys.argv[1] == "-t":
            file_format = FileFormat.TRNA
    else:
        data_file = open(sys.argv[1], 'r')

    if file_format == FileFormat.SS:
        trna_records = trnapy.ReadTRNARecordsFromSSFile(data_file)
    elif file_format == FileFormat.FA:
        trna_records = trnapy.ReadTRNARecordsFromFAFile(data_file)
    elif file_format == FileFormat.TRNA:
        trna_records = trnapy.ReadTRNARecordsFromTRNAFile(data_file)

    data_file.close()
    trna_records.AddCCAToAll()
    trna_records.ExpandAll()

    def PrintPatternRecord(trna):
        record_to_print = ""
        if trna.is_virtual:
            record_to_print += "!"
        #record_to_print += trna.sign
        if not trna.is_virtual:
            if trna.expanded:
                record_to_print += "1-"
            else:
                record_to_print += "0-"
            record_to_print += str(len(trna) - 4)
        else:
            record_to_print += "3-"
            if trna.expanded:
                record_to_print += str(len(trna) - 2)
            else:
                record_to_print += str(len(trna) - 1)

        record_to_print += " " + trna.sequence()
        print(record_to_print)

    for trna in trna_records:
        PrintPatternRecord(trna)
        PrintPatternRecord(trna.inverse_complement())

