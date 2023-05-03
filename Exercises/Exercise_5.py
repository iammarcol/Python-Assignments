# ITERATOR - used in all of the part5 exercises

def FASTA_iterator(fasta_filename):
    identif = ""
    sequence = ""
    with open(fasta_filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                if identif and sequence:
                    yield identif, sequence
                identif = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
    if identif and sequence:
        yield identif, sequence

###############################################################################

# PROBLEM_1a (session2_1)

def get_proteins_ratio_by_residue_threshold(fasta_filename,residue,relative_threshold=0.03,absolute_threshold=10):
    a=[]
    s=[]
    sequence=""
    abs_count=0
    rel_count=0
    prot_that_accomp=0
    tot_aa=0
    for identif, seq in FASTA_iterator(fasta_filename):
        a.append(identif)
        s.append(seq)
    num_prot=len(a)
    for sequence in s:
        tot_aa=len(sequence)       
        abs_count=sequence.count(residue)
        rel_count=float(abs_count/tot_aa)
        if abs_count>=absolute_threshold and rel_count>=relative_threshold:
            prot_that_accomp+=1
    ratio=float(prot_that_accomp/num_prot)
    print(ratio)   

# get_proteins_ratio_by_residue_threshold("example.fasta","A")

# PROBLEM_1b (session2_2)

def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
    
    ofd = open(output_filename, "w")
    # put the FASTA_iterator 
    for identifier, cseq in FASTA_iterator(filename):
        done_aminoacids = ""
        counts = ""
        for aminoacid in cseq:
            if aminoacid not in done_aminoacids:
                counts += "%s:%s," %(aminoacid, cseq.count(aminoacid))
                done_aminoacids += aminoacid
        ofd.write("%s\t%s\t%s\t%s\n" %(identifier[1:].strip(), 
                                       cseq[:first_n], 
                                       cseq[-last_m:], 
                                       counts.strip(",")))
    ofd.close()

# print_sequence_summary("example.fasta","ouuuuut.txt",first_n=10,last_m=10)


# PROBLEM_2

def get_max_sequence_length_from_FASTA_file(fasta_filename):
    return max([len(sequence[1]) for sequence in FASTA_iterator(fasta_filename)])


# PROBLEM_3

def get_min_sequence_length_from_FASTA_file(fasta_filename):
    return min([len(sequence[1]) for sequence in FASTA_iterator(fasta_filename)])


# PROBLEM_4

def get_longest_sequences_from_FASTA_file(fasta_filename):
    # put the output from the iterator into a single list
    fasta_list = list(FASTA_iterator(fasta_filename))   
    # find what is the length of the LONGEST seq.
    max_length = max(len(entry[1]) for entry in fasta_list)
    # put in a list if the length of seq. from fasta_list is = to the maximum length
    max_length_sequences = [entry for entry in fasta_list if len(entry[1]) == max_length]
    # sort
    max_length_sequences.sort(key=lambda x: x[0].lower())
    return max_length_sequences

# PROBLEM_5

def get_shortest_sequences_from_FASTA_file(fasta_filename):
    # put the output from the iterator into a single list
    fasta_list = list(FASTA_iterator(fasta_filename))   
    # find what is the length of the SHORTEST seq.
    min_length = min(len(entry[1]) for entry in fasta_list)
    # put in a list if the length of seq. from fasta_list is = to the minimum length
    min_length_sequences = [entry for entry in fasta_list if len(entry[1]) == min_length]
    # sort
    min_length_sequences.sort(key=lambda x: x[0].lower())
    return min_length_sequences

# PROBLEM_6
   
def get_molecular_weights(fasta_filename):
    dictionary={}
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    for identif, seq in FASTA_iterator(fasta_filename):
        result = 0
        for AA in seq:
            if AA in aminoacid_mw:
                result += aminoacid_mw[AA]
        dictionary[identif] = float(result)
    return dictionary

# PROBLEM_7

def get_sequence_with_min_molecular_weight(fasta_filename):
    # it uses previously defined funcion
    molecular_weights = get_molecular_weights(fasta_filename)
    # find the minimum value
    min_weight = min(molecular_weights.values()) 
    # compare it with the tuples generated from the iterator
    for identif, seq in FASTA_iterator(fasta_filename):
        if molecular_weights[identif] == min_weight:
            return (identif, seq)

# get_sequence_with_min_molecular_weight("input_file.fasta")

# PROBLEM_8

def get_mean_molecular_weight(fasta_filename):
    # it uses previously defined funcion
    molecular_weights = get_molecular_weights(fasta_filename)
    weights_list=[]
    total=0
    # find the minimum value
    for key,value in molecular_weights.items():
        # make a list with all values of the weights
        weights_list.append(value)
    for weight in weights_list:
        # sum all the values that are in the list
        total+=weight
    # use len of the weights_list, because that's the total num. of weights (seq)
    mean=total/len(weights_list)
    print(mean)

# get_mean_molecular_weight("input_file.fasta")