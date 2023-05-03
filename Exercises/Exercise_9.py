# PROBLEM_1-7_7(1-4)

import sys
import os

def errprint(text_string):
    sys.stderr.write(text_string)
    sys.stderr.flush()

###################################################################

class Sequence:
    alphabet = set()

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence.upper()

        if not set(self.__sequence).issubset(self.alphabet):
            invalid_chars = set(self.__sequence) - self.alphabet
            raise IncorrectSequenceLetter(invalid_chars,self.__class__.__name__)
            
    def get_identifier(self):
        return self.__identifier

    def get_sequence(self):
        return self.__sequence
        
    def get_mw(self):

        protein_weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        rna_weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
        dna_weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}

        if isinstance(self, ProteinSequence):
            weigths_dict = protein_weights
        elif isinstance(self, DNASequence):
            weigths_dict = dna_weights
        elif isinstance(self, RNASequence):
            weigths_dict = rna_weights
        
        return sum(weigths_dict[monomer] for monomer in self.get_sequence())

    def has_subsequence(self, sequence_obj):
        return sequence_obj.get_sequence() in self.__sequence

    def __len__(self):
        return len(self.get_sequence())

    def __eq__(self, other):
        return self.get_sequence() == other.get_sequence()

    def __ne__(self, other):
        return self.get_sequence() != other.get_sequence()

    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError("Can only combine same type of sequences!")
        
        comb_ident = self.get_identifier() + "+" + other.get_identifier()
        comb_seq = self.get_sequence() + other.get_sequence()

        return self.__class__(comb_ident, comb_seq)

    def __getitem__(self,key):
        return self.get_sequence()[key]

    def __contains__(self, item):
        return item in self.get_sequence()

    def __lt__(self, other):
        return self.get_mw() < other.get_mw() 
    def __le__(self, other):
        return self.get_mw() <= other.get_mw() 
    def __gt__(self, other):
        return self.get_mw() > other.get_mw() 
    def __ge__(self, other):
        return self.get_mw() >= other.get_mw() 

    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))


###################################################################

class NucleotideSequence(Sequence):
    
    def translate(self):

        rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
        rna_stop_codons = ['UAA', 'UAG', 'UGA']
        rna_start_codons = ['UUG', 'CUG', 'AUG']

        dna_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}
        dna_stop_codons = ['TAA', 'TAG', 'TGA']
        dna_start_codons = ['TTG', 'CTG', 'ATG']

        if isinstance(self, DNASequence):
            codon_table_dict = dna_table
            start_list = dna_start_codons
            stop_list = dna_stop_codons
        elif isinstance(self, RNASequence):
            codon_table_dict = rna_table
            start_list = rna_start_codons
            stop_list = rna_stop_codons
        else:
            raise ValueError('Invalid sequence type')

        sequence = self.get_sequence()
        aa_list = []

        has_start_index = [sequence.find(codon) for codon in start_list \
            if codon in sequence]
        
        if has_start_index:
            start_index = min([sequence.find(codon) for codon in start_list \
                if codon in sequence])
            for i in range(start_index, len(sequence), 3):
                codon = sequence[i: i+3]
                if codon in stop_list or len(codon) != 3 :
                    break
                aa_list.append(codon_table_dict[codon])
            
            return ''.join(aa_list)

        else:
            raise ValueError('Object has no start codon')

###################################################################

class DNASequence(NucleotideSequence):
    alphabet = set('GATC')
    
    def transcribe(self):
        return ''.join(['U' if D == 'A' else 'G' if D == 'C' else 'C' if D == 'G' else 'A' for D in self.get_sequence()])

    
###################################################################

class RNASequence(NucleotideSequence):
    alphabet = set('GAUC')
    
    def reverse_transcribe(self):
        return self.get_sequence().translate(str.maketrans('UGCA', 'ACGT'))

###################################################################

class ProteinSequence(Sequence):
    protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
    alphabet = set(protein_letters)

###################################################################

# creating a new ValueError
class IncorrectSequenceLetter(ValueError):
    def __init__(self, invalid_chars, class_name):
        self.invalid_chars = invalid_chars
        self.class_name = class_name
    
    def __str__(self):
        return f"The sequence item {self.invalid_chars} is not found in the alphabet of class {self.class_name}"

###################################################################

def FASTA_iterator(fasta_filename, SeqType):
    with open(fasta_filename) as file:
        sequence = ''
        # iterating through each line in the file
        for line in file:
            # check if the line contains fasta ID starting symbol
            if line.startswith('>'):
                if sequence:
                    try:
                        yield SeqType(identifier, sequence)
                    except IncorrectSequenceLetter as ISL:
                        errprint(f'{ISL}\n')
            # strip, concatinate
                identifier, sequence = line[1:].strip(), ''
            else:
                sequence += line.strip()
        if sequence:
            try:
        # if the creation of the seq obj. raises an "ISL" error
        # the exception is caught and an error message is printed
                yield SeqType(identifier, sequence)
            except IncorrectSequenceLetter as ISL:
            # prints the error message
                errprint(f'{ISL}\n')


###################################################################

def main():
    if len(sys.argv) not in [1, 2, 3]:
        print(f"Usage: \tpython3 {sys.argv[0]} [IN] [OUT] or\n"
              f"\tpython3 {sys.argv[0]} [IN] or\n"
              f"\tpython3 {sys.argv[0]}\n")
        sys.exit(1)
# output_path is set to "None" at first
# checking the number of comman-line arguments passed to the script
    input_path = sys.argv[1] if len(sys.argv) > 1 else '.'
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
# checking if input_path is a directory or a file
    files = []
    if os.path.isdir(input_path):
        files = [os.path.join(input_path, f) for f in os.listdir(input_path)
                 # search for files w the extension "fasta/fa"
                 if f.endswith(('.fasta', '.fa'))]
    # if it's a file it adds the full path 
    elif os.path.isfile(input_path):
        files = [input_path]

    if not files:
        print("Couldn't load any files")
        sys.exit(1)

    results = []
    for file in files:
        # using FASTA_iterator to loop through eachfile and reads the DNA seq.
        for seq in FASTA_iterator(file, DNASequence):
            # translating to protein seq.
            prot = ProteinSequence(seq.get_identifier(), seq.translate())
            # gets mw
            ident_len_mw = [prot.get_identifier(), len(prot), prot.get_mw()]
            # store everything in the "results" list
            results.append(ident_len_mw)

    print(f'{len(files)} FASTA files found.')
    print(f'{len(results)} sequences found.')
    print("Sorting the sequences...")
    # after finishing sort the results list by mw in ascending order
    sorted_results = sorted(results, key=lambda x: x[2], reverse=False)
    print("Sort process finished.")

    if output_path:
        with open(output_path, 'w') as f:
            for seq in sorted_results:
                f.write(f"{seq[0]}\t{seq[1]}\t{seq[2]}\n")
    else:
        for seq in sorted_results:
            print(f"{seq[0]}\t{seq[1]}\t{seq[2]}")
    # prints message that the program has finished correctly
    print("Program finished correctly.")
    # exit with a status code of 0 (successful)
    sys.exit(0)

if __name__ == '__main__':
    main()


