#!/usr/bin/python3
# Produces the *matured.fa file from the .ss file.
#
# The records of the *matured.fa file have the following format.
# <trna name>_<amino acid><anticodon>_<chromosome number>_<+->
#   _<start position>_<end position>

import re
import sys
import trnapy

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as ss_file:
        trna_records = trnapy.ReadTRNARecordsFromSSFile(ss_file)
    trna_records.AddCCAToAll()
    for trna in trna_records:
        chromosome_number = re.sub(r"chr([^.]*).*", r"\1",
                trna.identifier)
        trna_name = re.sub(r"chr[^.]*\.(.*)", r"\1", trna.identifier)
        matured_string = ">"
        matured_string += trna_name + "_"
        matured_string += trna.amino_acid + trna.anticodon + "_"
        matured_string += chromosome_number + "_"
        matured_string += trna.sign + "_"
        matured_string += str(trna.genome_interval[0]) + "_"
        matured_string += str(trna.genome_interval[1])
        print(matured_string)
        print(trna.sequence())

