# PROBLEM_1

### the idea
# make a list of all prot names, then len is the num_prot (the total num of proteins in the file)
# count the num of proteins that accomplish required input values for the thresholds
# extract sequences, calculate relative and absolute frequency

def get_proteins_ratio_by_residue_threshold(filename,residue,relative_threshold=0.03,absolute_threshold=10):
    a=[]
    s=[]
    sequence=""
    abs_count=0
    rel_count=0
    prot_that_accomp=0
    tot_aa=0
    with open(filename,"r") as ph:     # proteini=a    # sekvnence zjd=s
        for line in ph:
            if line.startswith('>'):
                a.append(line.replace("\n",""))
                if sequence!="":
                     s.append(sequence)
                     sequence=""
            else: 
                sequence+=line.replace("\n","")
        s.append(sequence)
        num_prot=len(a)
        for sequence in s:
            tot_aa=len(sequence)       
            abs_count=sequence.count(residue)
            rel_count=float(abs_count/tot_aa)
            if abs_count>=absolute_threshold and rel_count>=relative_threshold:
                prot_that_accomp+=1
        ratio=float(prot_that_accomp/num_prot)
        print(ratio)      


# PROBLEM_2

def print_sequence_summary(filename,output_filename,first_n=10,last_m=10):
    sequence=""
    ident=""
    line=""
    abs_count=0
    unique=[]  
    out=open(output_filename,'w')
    with open(filename,"r") as ph:
        for line in ph:
            if line.startswith('>'):
                string1="\n"
                out.write(string1)
                ident=line.replace("\n","")
                sequence=""  
            else:
                sequence+=line.replace("\n","")
                string2=f"{ident.lstrip('>')} \t {sequence[0:first_n]} \t {sequence[-last_m:]}"+"\t"
                out.write(string2)
                for aminoacid in sequence:     
                    abs_count=sequence.count(aminoacid)
                    if aminoacid not in unique:
                        unique.append(aminoacid)
                        string3=(f"{aminoacid}:{abs_count}"+',')
                        out.write(string3)
                unique=[]
    out.close()





