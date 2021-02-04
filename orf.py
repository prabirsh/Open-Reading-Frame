remove the new line character
def remove_newlines(fname):
    flist = open(fname).readlines()
    return [s.rstrip('\n') for s in flist]

line = remove_newlines("prabir_input.txt")

i = 0
a = str()
while i < len(line):
    a  += line[i]
    dna = a
    i += 1  

print ("length of the DNA:-",len(dna))
rna = ""
#generate RNA string
for i in dna:
    # Replace all occurrences of T with U
    if i == "T":
        rna += "U"
    else:
        rna += i

#generate rna file 
fout = open("prabir_output1_RNA.txt","w")
fout.writelines(rna)
fout.close()

#open the rna file and read the data
rnaopen = open("prabir_output1_RNA.txt","r")
print ("\n -::::::::::::::::RNA:::::::::::::::-")
print (">RNA Sequence from base 0 to ",len(dna))
for line in rnaopen:
    print (line)

print ("\n \n :::-File:prabir_output1_RNA.txt created on your directory where RNA sequences are stored-::: \n \n")

def codon(seq):
    #codon dict
    rna_codon = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
           "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
           "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
           "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
           "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
           "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }


    print ("\n-::::::::::::::Protein:::::::::::::::-")
    
    # Generate protein string
    print ("Enter 1 :-for Reading Frame 1 on the direct starnd")
    print ("Enter 2 :-for Reading Frame 2 on the direct starnd")
    print ("Enter 3 :-for Reading Frame 3 on the direct starnd")
    print ("Enter your Reading Frame:")
    reading_1 = int(input()) #choose reading frame
    if reading_1 == 1:
        print (">Protein Sequence from base 0 to ",len(dna),"bp")
        for i in range(0, len(line)-(3+len(line)%3), 3):
            codon = rna_codon[line[i:i+3]]
            protein(codon)
    elif reading_1 == 2:
        print (">Protein Sequence from base 0 to ",len(dna),"bp")
        for i in range(1, len(line)-(3+len(line)%3), 3):
            codon = rna_codon[line[i:i+3]]
            protein(codon)
    elif reading_1 == 3:
        print (">Protein Sequence from base 0 to ",len(dna),"bp")
        for i in range(2, len(line)-(3+len(line)%3), 3):
            codon = rna_codon[line[i:i+3]]
            pro_3 = protein(codon)
    else:
        print ("INPUT//Invalid!!")

def protein(codon):
    if codon == "STOP":
        print (end="*")
    else:
        print (codon,end = "")
        #protein output file
        proteinoutput = open("prabir_output2_protein.txt","a")
        proteinoutput.writelines(codon)
        proteinoutput.close()

codon = codon(line)         #function call
print ("\n \n :::-File:prabir_output2_protein.txt created on your directory where Protein sequences are stored-::: \n \n")
