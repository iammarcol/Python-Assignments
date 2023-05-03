# PROBLEM_part3

def calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions,output_filename):
    a=[]
    s=[]
    b=[]
    alone_seq=[] 
    prop=[]
    sequence=""
    seqdict={}
    sorted_dict={}
    out=open(output_filename,'w')
    with open(fasta_filename,"r") as ph:     # proteins=a    # sequences together=s
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
        for i in s:
            alone_seq.append(i)      # separated sequences = alone_seq
            alone_seq=[]
    with open(subsequences_filename,"r") as pt:       # subsequences
        for linee in pt:
            b.append(linee.replace("\n",""))
        num_subseq=len(b) 
        dist=len(max(b, key=len))  # find out what is the longest subseq for the alignment at the end
        for i in b:
            ii=number_of_repetitions*i
            for j in s:
                if ii in j:
                    prop.append(j)
            seqdict[i]=[len(prop),round((len(prop)/num_prot),4)]   # makes a dict
            prop=[]
        sorted_dict=sorted(seqdict.items(), key=lambda x: x[1][1],reverse=True)
        string1="#Number of proteins:"+f"{num_prot:>24}"+"\n"+"#Number of subsequences:"+f"{num_subseq:20}"+"\n"
        out.write(string1)
        out.write("#Subsequence proportions:"+"\n")
        for tuple_a in sorted_dict:
            primero=tuple_a[0]
            secondo=tuple_a[1]
            secondo_a=secondo[0]
            secondo_b=secondo[1]
            string2=f"{primero:<{dist}}\t{secondo_a:>12}\t{secondo_b:>20.4f}\n"
            out.write(string2)
    out.close()


# check: calculate_aminoacid_frequencies("input.fasta","subseq.fasta",3,"ex3out.txt")



