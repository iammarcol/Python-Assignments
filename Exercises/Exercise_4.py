# PROBLEM_1

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

# PROBLEM_2

# input to check
# fasta_filenames_list=["input.fasta","input2.fasta","input3.fasta"]

def compare_fasta_file_identifiers(fasta_filenames_list):
    a=[]
    fasta_id=[]
    filter=[]
    intersec=[]
    freq={}
    dict={}
    spec={}
    for i in fasta_filenames_list:
        it = FASTA_iterator(i) 
        for n in it:
            a.append(n[0].upper())     # make a list of all identifiers of one file and set to upper each ID for the case-insensivity
        fasta_id.append(a)     # make a list of all file's identifiers, for each file as a list (a)
        a=[]
    for each in fasta_id:
        for item in set(each):   # using 'set' if there's two of the same IDs in the same file
            filter.append(item) 
    for id in filter:
        freq[id]=filter.count(id)               # count how many files contain certain "id" 
        if filter.count(id) == len(fasta_id):   # if ID is in each file (represented as a list inside of a fasta_id list), then its count is the number of files 
            intersec.append(id)
    for i in fasta_filenames_list:
        spec[i]=set(fasta_id[fasta_filenames_list.index(i)]) - set(intersec)  # subtract from the intersection, because those can't be unique
        spec_list = list(spec[i])
        unique_ids = [id for id in spec_list if filter.count(id) == 1]   # filtering: if it counts 1 ID that means it's unique
        spec[i]=set(unique_ids)                                          # updates for the final dictionary based on the unique_ids for each file {i}
    dict={'intersection':set(intersec),'union':set(filter),'frequency':freq,'specific':spec} 
    print(dict)


# compare_fasta_file_identifiers(fasta_filenames_list)

###########################################################################

# I have also made the iterator this way, making a list with tuples
# and then iterating throught the list
# but I think the other iterator is more suitable for the exercise task

#def FASTA_iterator(fasta_filename):
#    a=[]
#    s=[]
#    sequence=""
#    with open(fasta_filename,"r") as ph:     # proteins=a  # sequences together=s
#        for line in ph:
#            if line.startswith('>'):
#                 a.append(line.replace("\n",""))
#                 if sequence!="":
#                     s.append(sequence)
#                     sequence=""
#            else: 
#                sequence+=line.replace("\n","")
#        s.append(sequence)
#        fasta=list(zip(a,s))
#        for tuple in fasta:
#            yield tuple

