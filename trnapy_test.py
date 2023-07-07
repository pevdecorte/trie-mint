#!/usr/bin/python3

import copy
import functools
import io

import testing
import trnapy

ut = testing.UnitTestCollection()

trna = trnapy.tRNARecord(
        "chr19.trna1",
        "GTCTCTGTGGCGCAATCGGTTAGCGCGTTCGGCTGTTAACCGAAAGGTTGGTGGTTCGAGCCCACCCAGGGACG",
        (1383562, 1383635), "+", "Asn", "GTT")

trna2 = trnapy.tRNARecord(
        "chr11.trna4",
        "GCCCGGATAGCTCAGTCGGTAGAGCATCAGACTTTTAATCTGAGGGTCCGGGGTTCAAGTCCCTGTTCGGGCG",
        (59323902, 59323974), "+", "Lys", "TTT")

trna3 = trnapy.tRNARecord(
        "chr8.trna11",
        "GTAGTCGTGGCCGAGTGGTTAAGGCGATGGACTAGAAATCCATTGGGGTCTCCCCGCGCAGGTTCGAATCCTGCCGACTACG",
        (96281885, 96281966), "-", "Ser", "AGA")

fragmentable_trna = trnapy.FragmentableTRNARecord(
        "chrX.trna1530", 33,
        "GGTCTCATGGTGTAATGGTtAGCACACTGGACTTTGAGTCCAGCAaTCCGAGTTCGAGTCTTGGTGAGACCA",
        (13859166, 13859237), "-", "Gln", "TTG")

trna_with_intron = trnapy.FragmentableTRNARecord(
    "chr7.trna259", 34,
    "GCTCCAGTGGCGCAATCGGTtAGCGCGCGGTACTTATAtgtcagtgctaggcgtaagcgATGCCGAGGtTGTGAGTTCGATCCTCACCTGGAGCA".upper(),
    (28372999, 28373093), "+", "Ile", "TAT", intron_interval=(39, 59))

@ut(trna)
def tRNARecords_add_cca_test(trna):
    ut.ExpectEq(trna.cca_added, False)
    ut.ExpectEq(trna.cca_added, False)
    no_cca_sequence = copy.copy(trna.sequence())
    trna.add_cca()
    ut.ExpectEq(trna.sequence(), no_cca_sequence + "CCA")
    ut.ExpectEq(trna.cca_added, True)

@ut(trna)
def trna_not_changed_test(trna):
    ut.AssertEq(trna.cca_added, False)

@ut(trna)
def tRNARecords_expand_test(trna):
    expanded_list = trna.expand()
    for t in expanded_list:
        ut.ExpectEq(t.identifier[:-2], trna.identifier)
        ut.ExpectEq(t.sequence()[1:], trna.sequence())
        ut.ExpectEq(t.expanded, True)
    ut.ExpectEq(trna.expanded, False)

@ut(trna)
def tRNARecords_fragment_test(trna):
    trna.add_cca()
    length = len(trna)
    fragment = trna.fragment(0, length)
    ut.ExpectEq(fragment.sequence(), trna.sequence())
    fragment = trna.fragment(5, 10)
    ut.ExpectEq(fragment.sequence(), "TGTGGCGCAA")
    ut.ExpectEq(fragment.source_trna, trna)

@ut(trna, trna2, trna3)
def tRNARecordsFrame_test(trna1, trna2, trna3):
    frame = trnapy.tRNARecordsFrame()
    frame.add(trna1)
    frame.add(trna2)
    frame.add(trna3)
    ut.ExpectEq(len(frame), 3)

    virtual_frame = frame.InverseComplements()
    ut.ExpectEq(len(virtual_frame), 3)
    for trna1, trna2 in zip(frame, virtual_frame):
        ut.ExpectEq(trna1.is_virtual, False)
        ut.ExpectEq(trna2.is_virtual, True)
        ut.ExpectEq(
                trna1.sequence(),
                trna2.inverse_complement().sequence())
        ut.ExpectEq(trna2.sign, "-" if trna1.sign == "+" else "+")

    frame.AddCCAToAll()
    frame.ExpandAll()

    ut.ExpectEq(len(frame), 15)



@ut(trna, trna2, trna3)
def tRNARecordsFrame_CCA_test(trna1, trna2, trna3):
    frame = trnapy.tRNARecordsFrame()
    frame.add(trna1)
    frame.add(trna2)
    frame.add(trna3)
    frame.AddCCAToAll()
    for trna in frame:
        ut.ExpectEq(trna.sequence()[-3:], "CCA")


@ut(trna, trna2, trna3)
def tRNARecordsFrame_expand_test(trna1, trna2, trna3):
    frame = trnapy.tRNARecordsFrame()
    frame.add(trna1)
    frame.add(trna2)
    frame.add(trna3)
    frame.AddCCAToAll()
    frame.ExpandAll()

    ut.ExpectEq(len(frame), 15)
    
    expanded_trnas = [ t for t in frame if t.expanded ]

    ut.ExpectEq(len(expanded_trnas), 12)

    for trna in frame:
        ut.ExpectEq(trna.sequence()[-3:], "CCA")


@ut(fragmentable_trna)
def FragmentableTRNARecord_test(trna):
    trna.add_cca()
    # Anticodon is at 33-35.
    fragment = trna.fragment(0, 35)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.FIVE_PRIME_TRH)
    fragment = trna.fragment(0, 36)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.FIVE_PRIME_TRF)
    fragment = trna.fragment(33, len(trna) - 33)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.THREE_PRIME_TRH)
    fragment = trna.fragment(33, len(trna) - 1 - 33)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.THREE_PRIME_TRH)
    fragment = trna.fragment(33, len(trna) - 2 - 33)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.THREE_PRIME_TRH)
    fragment = trna.fragment(33, len(trna) - 3 - 33)
    ut.ExpectEq(fragment.fragment_type(),
            trnapy.MatchType.I_TRF)

    frame = trnapy.tRNARecordsFrame()
    frame.add(trna)

    # Checks that "CCA" does not get added twice.
    sequence_length = len(frame["chrX.trna1530"])
    frame.AddCCAToAll()
    frame.ExpandAll()
    ut.ExpectEq(len(frame["chrX.trna1530"]), sequence_length)

    # Checks that anticodon indicies are moved over by 1
    # after expansion.
    ut.ExpectEq(frame["chrX.trna1530"].anticodon_start, 33)
    ut.ExpectEq(frame["chrX.trna1530_A"].anticodon_start, 34)


@ut()
def InverseComplement_test():
    ut.ExpectEq(trnapy.InverseComplement("ACTG"), "CAGT")


@ut()
def ReadTRNARecordsFromSSFile_test():
    string_from_file = """chr7.trna147 (19301244-19301330)	Length: 87 bp
Type: SeC	Anticodon: TCA at 36-38 (19301279-19301281)	Score: 146.4
HMM Sc=0.00	Sec struct Sc=0.00
         *    |    *    |    *    |    *    |    *    |    *    |    *    |    *    |    * 
Seq: GCCCGGATGATCCTCAGTGGTCTGGGGTGCAGGCTTCAAACCTGTAGCTGTTTAGCGACAGAGTGGTTCAATTCCACCTTTCGGGCG
Str: >>>>>>>.>..>>>>>>....<<<<<<>>>>>>.......<<<<<<.>>>>>....<<<<<.>>>>.......<<<<<.<<<<<<<.

chr7.trna259 (28372999-28373093)	Length: 95 bp
Type: Ile	Anticodon: TAT at 35-37 (28373033-28373035)	Score: 73.7
Possible intron: 39-59 (28373037-28373057)
HMM Sc=45.50	Sec struct Sc=28.20
         *    |    *    |    *    |    *    |    *    |    *    |    *    |    *    |    *    |    
Seq: GCTCCAGTGGCGCAATCGGTtAGCGCGCGGTACTTATAtgtcagtgctaggcgtaagcgATGCCGAGGtTGTGAGTTCGATCCTCACCTGGAGCA
Str: >>>>>>>..>>>>.........<<<<.>>>>>............................<<<<<.....>>>>>.......<<<<<<<<<<<<.

chr7.trna545 (58399315-58399386)	Length: 72 bp
Type: Glu	Anticodon: TTC at 34-36 (58399348-58399350)	Score: 71.2
HMM Sc=46.50	Sec struct Sc=24.70
         *    |    *    |    *    |    *    |    *    |    *    |    *    | 
Seq: TCCCACATGGTCTAGCGGTtAGGATTCCTGGTTTTCACCCAGGCGGCCCGGGTTCGACTCCCGGTGTGGGAA
Str: >>>>>>>..>>>>........<<<<.>>>>>.......<<<<<....>>>>>.......<<<<<<<<<<<<."""
    f = io.StringIO(string_from_file)
    frame = trnapy.ReadTRNARecordsFromSSFile(f)
    ut.ExpectEq(len(frame), 3)

    # Checks that introns are correctly deleted.
    ut.ExpectEq(frame["chr7.trna259"].sequence(),
        "GCTCCAGTGGCGCAATCGGTtAGCGCGCGGTACTTATAATGCCGAGGtTGTGAGTTCGATCCTCACCTGGAGCA".upper())

    # Checks that anticodon indices are changed from 1-based to
    # 0-based.
    ut.ExpectEq(frame["chr7.trna259"].anticodon_start, 34)


@ut()
def ReadTRNARecordsFromFAFile_test():
    string_from_file = """>trna17_HisGTG_M_+_12138_12206
GTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACC
>trna18_SerGCT_M_+_12207_12265
GAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCA
>trna19_LeuTAG_M_+_12266_12336
ACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTA
>trna20_GluTTC_M_-_14674_14742
GTTCTTGTAGTTGAAATACAACGATGGTTTTTCATATCATTGGTCGTGGTTGTAGTCCGTGCGAGAATA
>trna21_ThrTGT_M_+_15888_15953
GTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACA
>trna22_ProTGG_M_-_15956_16023
CAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGA"""
    frame = trnapy.ReadTRNARecordsFromFAFile(
            io.StringIO(string_from_file))
    ut.ExpectEq(len(frame), 6)
    ut.ExpectIn("chrM.trna21", frame)

@ut()
def MatchFromMatchFileData_test():
    string_from_ss_file = """chr6.trna258 (48205068-48205139)	Length: 72 bp
Type: Cys	Anticodon: GCA at 33-35 (48205100-48205102)	Score: 57.3
HMM Sc=38.70	Sec struct Sc=18.60
         *    |    *    |    *    |    *    |    *    |    *    |    *    | 
Seq: GGGGGTATAGCTCAGGGGTGGAGCATTTGACTGCAGATCAAGGGGtCCCTGTTTCAAATCCAGGTGCCCCCT
Str: >>>>>>>..>>>>.......<<<<.>>>>>.......<<<<<.....>>>>.........<<<<<<<<<<<."""

    frame = trnapy.ReadTRNARecordsFromSSFile(
            io.StringIO(string_from_ss_file))
    frame.AddCCAToAll()
    frame.ExpandAll()
    virtual_frame = frame.InverseComplements()

    string_from_match_file = """GCTCAGGGGTGGAGC 130581161 10 +chr6.trna258_G
GCTCAGGGGTGGAGC 130586323 9 +chr6.trna258
GCTCAGGGGTGGAGC 130586323 10 +chr6.trna258_A"""

    match_dict = {}
    trnapy.ReadMatchDictFromFile(io.StringIO(string_from_match_file),
            "chr10", frame, virtual_frame, match_dict)

    all_matches = functools.reduce(lambda x,y: x+y,
            match_dict.values())
    ut.ExpectEq(len(all_matches), 3)


@ut()
def ConstructTRNASpace_test():
    string_from_ss_file = """chr6.trna258 (48205068-48205139)	Length: 72 bp
Type: Cys	Anticodon: GCA at 33-35 (48205100-48205102)	Score: 57.3
HMM Sc=38.70	Sec struct Sc=18.60
         *    |    *    |    *    |    *    |    *    |    *    |    *    | 
Seq: GGGGGTATAGCTCAGGGGTGGAGCATTTGACTGCAGATCAAGGGGtCCCTGTTTCAAATCCAGGTGCCCCCT
Str: >>>>>>>..>>>>.......<<<<.>>>>>.......<<<<<.....>>>>.........<<<<<<<<<<<."""

    frame = trnapy.ReadTRNARecordsFromSSFile(
            io.StringIO(string_from_ss_file))

    trna_space = trnapy.ConstructTRNASpace(frame)

    ut.AssertIn(("chr6", "+"), trna_space)
    ut.AssertEq(len(trna_space[("chr6", "+")]), 1)
    ut.ExpectEq(trna_space[("chr6", "+")][0], (48205067, 48205138))
    

@ut()
def MatchFromMatchFileData_test2():
    string_from_fa_file = """>trna16_GlnTTG_12_-_50211190_50211264
TCTAGGATGTGGTGTGATAGGTAGCATGGAGAATTTTGGATTCTCAGGGGTAGGTTCAATTCCTATATTCTAGAA"""
    string_from_match_file = "GGAATTGAACCTACCCCTGAGAATCCA 20036620 15 !chr12.trna16_T"
    fields = string_from_match_file.strip().split()
    trna_records = trnapy.ReadTRNARecordsFromFAFile(
            io.StringIO(string_from_fa_file))
    trna_records.AddCCAToAll()
    trna_records.ExpandAll()
    virtual_trna_records = trna_records.InverseComplements()

    match = trnapy.MatchFromMatchFileData(
            "chr10", *fields, trna_records, virtual_trna_records)

    ut.ExpectEq(str(match), "chr12.trna16_T+!:15 chr10:20036620")


chr19_trna9_fa_file = """>trna9_LysTTT_19_-_41748142_41748214
GCCAGGATAGTTCAGGTGGTAGAGCATCAGACTTTTAACCTGAGGGTTCAGGGTTCAAGTCTCTGTTTGGGCG"""
chr19_trna9_match_file = "CATCAGACTTTTAACC 242907472 24 -chr19.trna9"
chr19_trna9_match_chr = "chr2"
chr19_trna9_genome_interval = (242907472, 242907487)

@ut(
    chr19_trna9_fa_file,
    chr19_trna9_match_file,
    chr19_trna9_match_chr,
    chr19_trna9_genome_interval)
def GenomePosition_test(
        string_from_fa_file,
        string_from_match_file,
        match_chromosome,
        expected_genome_interval):
    trna_records = trnapy.ReadTRNARecordsFromFAFile(
            io.StringIO(string_from_fa_file))
    trna_records.AddCCAToAll()
    trna_records.ExpandAll()
    virtual_trna_records = trna_records.InverseComplements()

    fields = string_from_match_file.strip().split()

    match = trnapy.MatchFromMatchFileData(
            match_chromosome, *fields, trna_records, virtual_trna_records)

    ut.ExpectEq(
            match.unmatured_genome_interval(fields[0]),
            expected_genome_interval)


chr16_trna18_fa_file = """>trna18_GlyGCC_16_+_70822597_70822667
GCATTGGTGGTTCAGTGGTAGAATTCTCGCCTGCCATGCGGGCGGCCGGGCTTCGATTCCTGGCCAATGCA
>trna13_GlyCCC_17_+_19764175_19764245
GCATTGGTGGTTCAATGGTAGAATTCTCGCCTCCCACGCAGGAGACCCAGGTTCGATTCCTGGCCAATGCA"""
chr16_trna18_match_file = "GATTCCTGGCCAATGCACC 209262122 54 +chr16.trna18"
chr16_trna18_match_chr = "chr2"
chr16_trna18_genome_interval = (209262122, 209262138)
@ut(
    chr16_trna18_fa_file,
    chr16_trna18_match_file,
    chr16_trna18_match_chr,
    chr16_trna18_genome_interval)
def GenomePosition_test2(
        string_from_fa_file,
        string_from_match_file,
        match_chromosome,
        expected_genome_interval):
    trna_records = trnapy.ReadTRNARecordsFromFAFile(
            io.StringIO(string_from_fa_file))
    trna_records.AddCCAToAll()
    trna_records.ExpandAll()
    virtual_trna_records = trna_records.InverseComplements()

    fields = string_from_match_file.strip().split()

    match = trnapy.MatchFromMatchFileData(
            match_chromosome, *fields, trna_records, virtual_trna_records)

    ut.ExpectEq(
            match.unmatured_genome_interval(fields[0]),
            expected_genome_interval)


chr2_trna20_fa_file = """>trna20_GluTTC_2_-_131094701_131094772
TCCCATATGGTCTAGCGGTTAGGATTCCTGGTTTTCACCCAGGTGGCCCGGGTTCGACTCCCGGTATGGGAA
>trna5_GluTTC_13_-_41634874_41634945
TCCCATATGGTCTAGCGGTTAGGATTCCTGGTTTTCACCCAGGTGGCCCGGGTTCGACTCCCGGTATGGGAA"""
chr2_trna20_match_file = "GTTCCCATACCGGGAGTCGAACCCGGGCCACCTGGGTGAAAA 41634872 2 !chr13.trna5_A"
chr2_trna20_match_chr = "chr13"
chr2_trna20_genome_interval = (41634873, 41634913)

@ut(
    chr2_trna20_fa_file,
    chr2_trna20_match_file,
    chr2_trna20_match_chr,
    chr2_trna20_genome_interval)
def GenomePosition_test2(
        string_from_fa_file,
        string_from_match_file,
        match_chromosome,
        expected_genome_interval):
    trna_records = trnapy.ReadTRNARecordsFromFAFile(
            io.StringIO(string_from_fa_file))
    trna_records.AddCCAToAll()
    trna_records.ExpandAll()
    virtual_trna_records = trna_records.InverseComplements()

    fields = string_from_match_file.strip().split()

    match = trnapy.MatchFromMatchFileData(
            match_chromosome, *fields, trna_records, virtual_trna_records)

    ut.ExpectEq(
            match.unmatured_genome_interval(fields[0]),
            expected_genome_interval)



if __name__ == "__main__":
    ut.RunTests()
