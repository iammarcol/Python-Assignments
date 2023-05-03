# PROBLEM_1_(1-8)

class Sequence:
# class attribute 
    alphabet=""
    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence
    def get_identifier(self):
        return str(self.__identifier)
    def get_sequence(self):
        return str(self.__sequence)  
    def get_mw(self):
    # checks for the correct instance and then calculates mw depending on that
        if isinstance(self, ProteinSequence):
            weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        elif isinstance(self, RNASequence):
            weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
        elif isinstance(self, DNASequence):
            weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
    # since it checks for the correct instance here too I added ValueError again
        return sum(weights[monomer] for monomer in self.get_sequence())
    def has_subsequence(self, sequence_obj):
        return sequence_obj.get_sequence() in self.__sequence
    def __len__(self):
        return len(self.__sequence)
    def __eq__(self, other):             # this has to be checked as print(prot1==prot2)
        if isinstance(other, Sequence):
            return self.__sequence == other.get_sequence()
        return False
    def __ne__(self, other):
        if isinstance(other, Sequence):
            return self.__sequence != other.get_sequence()
        return True
    def __getitem__(self, i):           # has to be checked as: print(dna.get_sequence()[4])
        return self.__sequence[i]
    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError("Operands must be of the same class.")
        identifier = self.get_identifier() + "+" + other.get_identifier()
        sequence = self.get_sequence() + other.get_sequence()
        return type(self)(identifier, sequence)
    def __contains__(self, substring):    # has to be checked as: print("ATC" in dna.get_sequence()) 
        return substring in self.__sequence  # returns a boolean
    # comparing sequences, if the list is made and sorted
    # the sequences will be sorted based on mw starting from the lowest mw
    def __lt__(self, other):
        return self.get_mw() < other.get_mw()
    def __le__(self, other):
        return self.get_mw() <= other.get_mw()
    def __gt__(self, other):
        return self.get_mw() > other.get_mw()
    def __ge__(self, other):
        return self.get_mw() >= other.get_mw()
    # adapt the sequence class so that it can be used as key in a dictionary or it can be added to a set
    def __hash__(self):
        return hash((self.__identifier, self.__sequence))
    
#########################################

class NucleotideSequence(Sequence):
    # from RNA to protein
    rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)
    def translate(self):
        protein_seq = ""
        alphabet='ACDEFGHIKLMNPQRSTVWY'
        for i in range(0, len(self.get_sequence()), 3):
            codon = self.get_sequence()[i:i+3].upper()
            if len(codon) < 3:
                continue
            if codon in NucleotideSequence.rna_table:
                amino_acid = NucleotideSequence.rna_table[codon]
            else:
                amino_acid = 'X'
            if amino_acid in alphabet:
                protein_seq += amino_acid
            else:
                protein_seq += 'X'
        return protein_seq
    
#########################################

class DNASequence(NucleotideSequence):
    alphabet='GATC'
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    def __init__(self, identifier, sequence):
        # Check that all characters in sequence are valid
        invalid_chars = set(sequence) - set(self.alphabet)
        if invalid_chars:
            raise ValueError(f"Impossible to create instance: {', '.join(invalid_chars)} not possible")
        super().__init__(identifier, sequence) # sets the sequence attribute
    def transcribe(self):
        rna_sequence = ''
        for nucleotide in self.get_sequence():  # use method to access private variable
            rna_sequence += DNASequence.dna_complement.get(nucleotide, 'N')
        return rna_sequence

#########################################

class RNASequence(NucleotideSequence):
    alphabet='GAUC'
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}
    def __init__(self, identifier, sequence):
        # Check that all characters in sequence are valid
        invalid_chars = set(sequence) - set(self.alphabet)
        if invalid_chars:
            raise ValueError(f"Impossible to create instance: {', '.join(invalid_chars)} not possible")
        super().__init__(identifier, sequence)
    def reverse_transcribe(self):
        dna_seq = ''
        for base in self.get_sequence():
            dna_seq += RNASequence.dna_complement.get(base, 'N')
        return dna_seq
    
#########################################

class ProteinSequence(Sequence):
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    def __init__(self, identifier, sequence):
        # Check that all characters in sequence are valid
        invalid_chars = set(sequence) - set(self.alphabet)
        if invalid_chars:
            raise ValueError(f"Impossible to create instance: {', '.join(invalid_chars)} not possible")
        super().__init__(identifier, sequence)