# This module provides some basic functionality for working
# with tRNA's.

import copy
import re
from enum import Enum, auto

class FileFormat(Enum):
    SS = 1
    FA = 2
    TRNA = 3


# Represents a tRNA.
#
# original_trna: Points back to the original tRNA if this tRNA
#   is virtual. Otherwise should be set to None.
class tRNARecord:
    def __init__(self, identifier, sequence,
            genome_interval, sign, amino_acid, anticodon,
            is_virtual=False, original_trna=None):
        self.identifier = identifier
        self._sequence = sequence.upper()
        self.genome_interval = genome_interval
        # sign indicates the strand from which the tRNA originated.
        self.sign = sign
        self.amino_acid = amino_acid
        self.anticodon = anticodon

        # This will be set to True if this is the inverse
        # complement of an original tRNA, not a tRNA itself.
        self.is_virtual = is_virtual

        # True if "CCA" has already been appended.
        self.cca_added = False  

        # True if this tRNA comes from another by prepending
        # A, C, T, or G.
        self.expanded = False

    def __str__(self):
        return self.identifier

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._sequence)

    def sequence(self):
        return self._sequence

    # Returns a list of 4 newly created tRNA's,
    # obtained by prepending A, C, T, or G to the original.
    def expand(self):
        # Ensures we don't try to expand an already expanded tRNA.
        assert not self.expanded
        expanded_records = []
        for c in ['A', 'C', 'T', 'G']:
            new_trna = copy.deepcopy(self)
            new_trna.identifier += "_%s" % c
            new_trna._sequence = c + new_trna._sequence
            new_trna.expanded = True
            expanded_records.append(new_trna)
        return expanded_records
    
    # Returns a fragment of this tRNA. No type is computed.
    def fragment(self, start_index, length):
        # Nonsense to make a fragment unless CCA has been added.
        assert self.cca_added
        return Fragment(self, start_index, length)

    # Appends "CCA" to the sequence and updates the cca_added bit.
    # Does nothing if the cca_added bit is already set.
    def add_cca(self):
        if self.cca_added:
            return
        self._sequence += "CCA"
        self.cca_added = True

    # Creates and returns a new virtual tRNA which is the
    # inverse complement of this one.
    #
    # Virtual tRNA's are not fragmentable, which is why this
    # function returns an ordinary tRNARecord.
    def inverse_complement(self):
        virtual_trna =tRNARecord(self.identifier,
                InverseComplement(self._sequence),
                self.genome_interval,
                "-" if self.sign == "+" else "+",
                self.amino_acid,
                self.anticodon,
                is_virtual = True,
                original_trna = self)
        virtual_trna.expanded = self.expanded
        virtual_trna.cca_added = self.cca_added
        return virtual_trna


# Holds several tRNARecord objects and allows for some operations,
# such as appending "CCA", to be performed on all of them at once.
# tRNA's can be looked up by their identifiers, which is required
# when creating Match objects.
class tRNARecordsFrame:
    def __init__(self):
        self.trna_records = {}

    def __len__(self):
        return len(self.trna_records)

    def __str__(self):
        return "<tRNARecordsFrame with " + str(len(self)) + " records>"

    def __repr__(self):
        return str(self)

    # Implements lookups by identifier.
    def __getitem__(self, identifier):
        if identifier not in self.trna_records:
            raise KeyError(identifier)
        return self.trna_records[identifier]

    def __contains__(self, identifier):
        return identifier in self.trna_records

    def __iter__(self):
        return iter(self.trna_records.values())

    def add(self, trna):
        if not isinstance(trna, tRNARecord):
            raise TypeError("not a tRNARecord object")
        assert trna.identifier not in self.trna_records
        self.trna_records[trna.identifier] = trna

    def AddCCAToAll(self):
        for trna in self.trna_records.values():
            trna.add_cca()

    def ExpandAll(self):
        new_trna_records = {}
        for identifier, trna in self.trna_records.items():
            new_trna_records[identifier] = trna
            for expanded_trna in trna.expand():
                new_trna_records[expanded_trna.identifier] = \
                        expanded_trna
        self.trna_records = new_trna_records

    def InverseComplements(self):
        virtual_frame = tRNARecordsFrame()
        for trna in self:
            virtual_frame.add(trna.inverse_complement())
        return virtual_frame



# Represents a tRNA as read from a .ss file.
#
# Fragmentable tRNA records are meant to be those tRNA recods
# whose fragments have "types", e.g. 5'-tRF, 3'-tRH.
# Fragments can also be taken of ordinary tRNA's, but they will not
# have types.
#
# The intronic subsequence is automatically deleted.
# The original sequence is stored in the original_sequence
# member variable.
#
# Fragments can be created by calling the fragment() method.
# The extra data carried by a tRNARecordFragmentable object
# allows the fragment type (e.g. 3'-tRH, i-tRF) to be computed.
#
# By design, virtual tRNAs (inverse complements) are not fragmentable.
#
# The anticodon_start index, and well as the indices in intron_interval,
# should be 0-based. (They will have
# to be converted from the 1-based indices in the .ss file.)
class FragmentableTRNARecord(tRNARecord):
    def __init__(self, identifier, anticodon_start,
            sequence, genome_interval, sign,
            amino_acid, anticodon, intron_interval = []):
        super().__init__(identifier, sequence, genome_interval,
                sign, amino_acid, anticodon)

        self.anticodon_start = anticodon_start
        # The anticodon region always has length 3.
        self.anticodon_end = anticodon_start + 2
        self.original_sequence = self._sequence

        # Deletes the intronic subsequence if there is one. 
        if intron_interval:
            self._sequence = \
                    self.original_sequence[:intron_interval[0]] +\
                    self.original_sequence[intron_interval[1] + 1:]
            if intron_interval[0] < anticodon_start:
                anticodon_start -=\
                        intron_interval[1] - intron_interval[0] + 1

    # Calls the same method from super(), but also updates
    # the position of the anticodon and the intron interval.
    def expand(self):
        expanded_records = super().expand()
        for trna in expanded_records:
            trna.anticodon_start += 1
            trna.anticodon_end += 1
        return expanded_records

    # Names for fragments are computed as follows.
    # 5'-tRH: position 0 (-1) to within anticodon region
    # and -1 or -2 position
    # of anticodon (from .ss file) [if A1A2A3 is the anticodon
    # triplet and n1n2A1A2A3n3n4 denotes the anticodon loop and
    # the sequence immediately
    # surrounding the triplet, then 5′-tRHs terminate at any of the four
    # underlined positions n1∇n2∇A1∇A2∇A3n3n4 (each ∇ denotes end)]
    # 0 or -1 to n1, n2, A1, or A2
    #
    # 3'-tRH: n2, A1, A2 or A3 to any position in CCA 
    #
    # 5'-tRFs: position 0 to any position except n1, n2, A1, or A2
    #
    # 3'-tRFs: any position except n2, A1, A2 or A3 to any position in CCA
    #
    # i-tRFs: start after position 0 and end before the first C of CCA
    # addition
    def FragmentType(self, start_index, length):
        assert 0 <= start_index < len(self)

        # Computes the type of this fragment. 
        end_index = start_index + length - 1

        assert 0 <= end_index < len(self)

        if start_index == 0:
            # 5'-tRH
            if self.anticodon_start - 2 <= end_index <=\
                   self.anticodon_end - 1:
                fragment_type = MatchType.FIVE_PRIME_TRH
            # 5'-tRF
            else:
                fragment_type = MatchType.FIVE_PRIME_TRF
        elif end_index >= len(self) - 3:  # Match ends inside CCA region.
            # 3'-tRH
            if self.anticodon_start - 1 <= start_index \
                    <= self.anticodon_end:
                fragment_type = MatchType.THREE_PRIME_TRH
            # 3'-tRF
            else:
                fragment_type = MatchType.THREE_PRIME_TRF
        # i-tRF
        elif start_index > 0 and end_index < len(self) - 3:
            fragment_type = MatchType.I_TRF
        else:
            fragment_type = None  # Something has gone wrong.

        return fragment_type
    
    def fragment(self, start_index, length):
        # Nonsense to make a fragment unless CCA has been added.
        assert self.cca_added
        return Fragment(self, start_index, length,
                fragment_type=self.FragmentType(start_index, length))


# Computes and returns the inverse complement of a fragment,
# represented as a string.
def InverseComplement(fragment_str):
    # Returns the complement of a nucleotide.
    def NucleotideComplement(c):
        if c == 'A':
            return 'T'
        elif c == 'T':
            return 'A'
        elif c == 'C':
            return 'G'
        elif c == 'G':
            return 'C'
        else:
            return c
    fragment_list = list(fragment_str)
    fragment_list.reverse()
    ic = []
    for c in fragment_list:
        ic.append(NucleotideComplement(c))
    return "".join(ic)


# Reads a .ss file and returns the all the records as a tRNA
# records frame. The values are tRNARecord objects.
#
# Positions in the .ss file are 1-based, which is why there are
# lots of -1's in the code below.
def ReadTRNARecordsFromSSFile(ss_file):
    trna_records = tRNARecordsFrame()
    try:
        while True:
            intron_present = False
            line = next(ss_file)
            fields = line.split()
            identifier = fields[0]
            
            genome_interval_string = fields[1].strip(")(")
            genome_interval = genome_interval_string.split('-')
            genome_interval = (int(genome_interval[0]) - 1,\
                    int(genome_interval[1]) - 1)
            if genome_interval[0] > genome_interval[1]:
                sign = "-"
                genome_interval = \
                        (genome_interval[1], genome_interval[0])
            else:
                sign = "+"

            line = next(ss_file)
            fields = line.strip().split()
            amino_acid = fields[1]
            anticodon = fields[3]

            # Anticodon indices are 1-based in the .ss file, but
            # the FragmentableTRNARecord constructor expects
            # 0-based indices.
            anticodon_interval = fields[5]
            anticodon_start = int(anticodon_interval.split('-')[0]) - 1
            line = next(ss_file)
            
            # Checks whether a possible intron is present.
            # Intron intervals are 1-based in the .ss file, but
            # the FragmentableTRNARecord constructor expects
            # 0-based indices.
            if line.split()[1] == "intron:":
                intron_interval = line.split()[2].split('-')
                intron_interval = [int(x) - 1 for x in intron_interval]
                intron_present = True
                line = next(ss_file)

            # Reads in the tRNA sequence.
            while line.split()[0] != "Seq:":
                line = next(ss_file)
            sequence = line.strip().split()[1]

            if intron_present:
                trna_records.add(
                    FragmentableTRNARecord(identifier, anticodon_start,
                        sequence, genome_interval, sign,
                        amino_acid, anticodon, intron_interval))
            else:
                trna_records.add(
                    FragmentableTRNARecord(identifier, anticodon_start,
                        sequence, genome_interval, sign,
                        amino_acid, anticodon))

            # Scans until the next blank line.
            while line.strip() != "":
                line = next(ss_file)

    except(StopIteration):
        pass

    return trna_records


# Reads in tRNA records from a .fa file.
# In the .fa files we've seen, the left and right endpoints
# of the tRNA position in the genome are not reversed in
# order to indicate that they are on the minus strand.
#
# Positions are 0-based in the .fa file, unlike in the .ss file.
def ReadTRNARecordsFromFAFile(fa_file):
    trna_records = tRNARecordsFrame()
    for line in fa_file:
        t_identifier, amino_acid_and_anticodon, \
            chromosome, sign, genome_interval_l, genome_interval_r = \
                    line.strip().split('_')
        sequence = next(fa_file).strip()
        t_identifier = t_identifier[1:]  # Removes leading '>'.
        identifier = "chr" + chromosome + "." + t_identifier
        amino_acid = amino_acid_and_anticodon[:3]
        anticodon = amino_acid_and_anticodon[3:]
        genome_interval_l = int(genome_interval_l) - 1
        genome_interval_r = int(genome_interval_r) - 1

        trna_records.add(tRNARecord(identifier, sequence,
            [genome_interval_l, genome_interval_r], sign,
            amino_acid, anticodon))

    return trna_records


# Reads in tRNA records from a *.trna file. This file format was created
# specifically for the trie-mint pipeline. The records in the *.trna file
# have the following format.
#
#   <trna descriptor> <amino acid> <anticodon> [+-] <start position> <end position>
#       <trna sequence>
#
# The start and end positions should be 0-based.
#
# For example:
#
#   chrM.trna18 Ser GCT + 12206 12264
#   GAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCA
def ReadTRNARecordsFromTRNAFile(trna_file):
    trna_records = tRNARecordsFrame()
    for line in trna_file:
        descriptor, amino_acid, anticodon, \
            sign, genome_interval_l, genome_interval_r = \
                    line.strip().split(' ')
        sequence = next(trna_file).strip()

        trna_records.add(tRNARecord(descriptor, sequence,
            [int(genome_interval_l) - 1, int(genome_interval_r) - 1], sign,
            amino_acid, anticodon))

    return trna_records

# Represents a tRNA fragment. The type of the fragment is
# automatically computed. The source trna is stored.
# The fragment_type field can be omitted, in which case
# it is set to None.
class Fragment:
    def __init__(self, source_trna, start_index, length,
            fragment_type=None):
        self.source_trna = source_trna
        self.start_index = start_index
        self.length = length

    def __str__(self):
        return self.sequence()

    def __repr__(self):
        return "<fragment from {}:{}-{} {}>".format(
                self.source_trna.identifier,
                self.start_index,
                self.start_index + self.length - 1,
                str(self))

    def __len__(self):
        return self.length

    def sequence(self):
        seq = self.source_trna.sequence()[
                self.start_index : self.start_index + self.length]
        if self.source_trna.is_virtual:
            return InverseComplement(seq)
        else:
            return seq

    # Only makes sense if source_trna is fragmentable.
    def fragment_type(self):
        assert isinstance(self.source_trna, FragmentableTRNARecord)
        return self.source_trna.FragmentType(
                self.start_index,
                self.length)


# Represents a partial match between a fragment (pattern) and part of the
# haystack.
# chromosome denotes the chromosome on which the match was found,
# not the chromosome from which the tRNA fragment originated.
#
# Prints a "!" after the identifier of a virtual tRNA.
#
# chromosome, strand_sign: Indicates where the match was found.
class Match(Fragment):
    def __init__(self, chromosome, strand_sign, haystack_start,
            source_trna, start_index, fragment_length):
        self.haystack_start = haystack_start
        self.chromosome = chromosome
        self.strand_sign = strand_sign
        super().__init__(source_trna, start_index, fragment_length)

    def __str__(self):
        return (
                self.source_trna.identifier + self.source_trna.sign +
                ("!" if self.source_trna.is_virtual else "") +
                ":" + str(self.start_index) + " " + self.chromosome +
                ":" + str(self.haystack_start))

    def __repr__(self):
        return str(self)

    # Returns an ordered pair (start, end) of genome indices
    # of the matched trna fragment, stripped of "CCA" and the
    # extra initial [ACTG]
    def unmatured_genome_interval(self, fragment_str=None):
        start_index = self.start_index
        haystack_start = self.haystack_start
        haystack_end = haystack_start + len(self) - 1
        end_index = start_index + len(self) - 1

        if self.source_trna.is_virtual: 
            if self.source_trna.expanded:
                if end_index == len(self.source_trna) - 1:
                    haystack_end -= 1
            if self.source_trna.cca_added:
                if start_index <= 2:
                    haystack_start += 3 - start_index
        else:
            if self.source_trna.expanded:
                if start_index == 0:
                    haystack_start += 1
            if self.source_trna.cca_added:
                if end_index >= len(self.source_trna) - 3:
                    haystack_end -= 4 - len(self.source_trna) + end_index

        return (haystack_start, haystack_end)



# Enumeration which represents the different types of matches.
class MatchType(Enum):
    FIVE_PRIME_TRH = auto()
    THREE_PRIME_TRH = auto()
    FIVE_PRIME_TRF = auto()
    THREE_PRIME_TRF = auto()
    I_TRF = auto()

    def __str__(self):
        if self.value == self.I_TRF.value:
            return "i-tRF"
        elif self.value == self.FIVE_PRIME_TRH.value:
            return "5'-tRH"
        elif self.value == self.THREE_PRIME_TRH.value:
            return "3'-tRH"
        elif self.value == self.FIVE_PRIME_TRF.value:
            return "5'-tRF"
        elif self.value == self.THREE_PRIME_TRF.value:
            return "3'-tRF"


def MatchFromMatchFileData(chromosome, fragment_str,
        haystack_start, start_index, source_trna_complex_name,
        trna_records, virtual_trna_records):
    haystack_start = int(haystack_start)
    start_index = int(start_index)

    stripped_trna_name = StripSignFromTRNAName(source_trna_complex_name)
    sign = ExtractSign(source_trna_complex_name)
    if sign == "!":
        source_trna = virtual_trna_records[stripped_trna_name]
    else:
        source_trna = trna_records[stripped_trna_name]

    return Match(
            chromosome, "-" if sign == "!" else "+",
            haystack_start, source_trna, start_index,
            len(fragment_str))


# Reads in matches from matches_file, and adds to match_dict.
# match_dict is a dictionary whose keys have the form
# (<chromosome>, <sign>, <(interval_start, interval_end)>).
# E.g. (chr19, "-", (160766110, 160766182)).
# 
# This tuple describes where the match was found, not the origin
# of the tRNA.
# Each value is a list of Match objects describing the matches
# at the given chromosome, strand (+ or -), and position.
def ReadMatchDictFromFile(matches_file, chromosome,
        trna_records, virtual_trna_records, match_dict):
    for line in matches_file:
        fields = line.strip().split(' ')
        match = MatchFromMatchFileData(
                chromosome, *fields,
                trna_records, virtual_trna_records)

        haystack_start = match.haystack_start
        haystack_end = match.haystack_start + len(match) - 1

        genome_interval = match.unmatured_genome_interval()
        key = (match.chromosome, match.strand_sign, genome_interval)

        if key not in match_dict:
            match_dict[key] = []
        match_dict[key].append(match)


# Returns a dictionary keyed by a <chromosome name, sign> pair.
# E.g. (chr19, "-")
# Each value is a list of ordered pairs <start_index, end_index>,
# describing the intervals making up tRNA space.
# The lists are sorted.
# The intervals on each strand are assumed to be pairwise disjoint.
def ConstructTRNASpace(trna_records):
    trna_space = {}
    for trna in trna_records:
        key = (ExtractChromosome(trna.identifier), trna.sign)
        if key not in trna_space:
            #trna_space[key] = []
            trna_space[key] = set()
        #trna_space[key].append(trna.genome_interval)
        trna_space[key].add(tuple(trna.genome_interval))

    for key, value in trna_space.items():
        trna_space[key] = list(value)

    for intervals in trna_space.values():
        intervals.sort(key=lambda x: x[0])

    return trna_space


######### A few functions for parsing trna names. #######
#
# trna_complex_name should have the form [+-]<trna_name>[_[ACTG]].
# Example: +chr1.trna702_G

# Returns: chr1.trna702_G
def StripSignFromTRNAName(trna_complex_name):
    m = re.search("[^-+!]+", trna_complex_name)
    return m.group(0)

# Returns: chr1
def ExtractChromosome(trna_complex_name):
    m = re.search("[^-+.]+", trna_complex_name)
    return m.group(0)

# Returns: +
def ExtractSign(trna_complex_name):
    return trna_complex_name[0]

# Returns: +chr1.trna702 
def DropExpansion(trna_complex_name):
    m = re.search("[^_]+", trna_complex_name)
    return m.group(0)

##########################################################
