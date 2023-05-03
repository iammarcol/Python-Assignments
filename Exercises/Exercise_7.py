# PROBLEM_1

class Protein:
    def __init__(self,identifier,sequence):
        self.identifier=identifier
        self.sequence=sequence
    def get_identifier(self):
        return str(self.identifier)
    def get_sequence(self):
        return str(self.sequence)
    def get_mw(self):
        aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        mw=0
        for AA in self.sequence:
            if AA in aminoacid_mw:
                mw+=aminoacid_mw[AA]
        return mw
    def has_subsequence(self,Protein):
        return Protein.get_sequence() in self.sequence
    def get_length(self):
        return len(self.sequence)

# PROBLEM_2

def FASTA_iterator(fasta_filename):
    identif = ""
    sequence = ""
    with open(fasta_filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                if identif and sequence:
                    yield Protein(identif, sequence)
                identif = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
    if identif and sequence:
        yield Protein(identif, sequence)



# check:

# if __name__=="__main__":
#     for protein in FASTA_iterator("test.fasta"):
#         print(protein.get_identifier())
#         print(protein.get_sequence())
#         print(protein.get_mw())
#         print(protein.has_subsequence(Protein("testghhhh","SATVSEINSETDFV")))
#         print(protein.get_length())
#         print()
