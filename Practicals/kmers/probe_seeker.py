#!/usr/bin/env python
# coding: utf-8

class sequence():
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

#Function parses file in .fasta format and returns a list of tuples in a form: 
#[(gene_name1, sequence1), (gene_name2, sequence2), ..., (gene_nameN, sequenceN)]
def parser():
    file = open("yeast.fa")
    
    sequence_lst = []
    full_sequence = ""
    sequence_name = ""
    
    for line in file:
        if line[0] == ">" and len(full_sequence) == 0:
            sequence_name = line[1::].split()[0]
        elif line[0] != ">":
            full_sequence += line.strip("\n")
        elif line[0] == ">" and len(full_sequence) > 0:
            sequence_lst.append(sequence(sequence_name,full_sequence))
            full_sequence = ""
            sequence_name = ""
            sequence_name = line[1::].split()[0]
            
    sequence_lst.append(sequence(sequence_name,full_sequence))
    
    file.close()
    
    return(sequence_lst)     

#Function returns k-mer spectrum of a given sequence (s) and given length (k)
def kmer_generator(s,k):
    lst = []
    for i in range(0,len(s) - k + 1):
        lst.append(s[i:i+k])
    return(lst)

long_seq_lst = parser()


start = 0
stop = 600


#The whole premise of the loop below is to create two sets of spectra. First one is 'global' spectrum 
#for all the given genes and then in each iteration of the given gene create a spectrum for all the gene 
#except this gene. After that generate an outcome of the difference beetween 'global spectrum' set and
#'global spectrum but one gene' set. This operation should give a set of unique sequences that only exist for
#this specific gene. By iterating for all the sequences for the given k will give us sets of such unique values
#for all the given sequences. If during Iterationthere is an empty set, increase the value of k-parametry by 1
#and exit the loop.

k = 3
dic = {}
counter = 0
seq_lst = long_seq_lst[start:stop]

while k < 100:
    print()
    print("##########################")
    print("probe length: ", k)
    all_kmers = set([kmer for seq in seq_lst for kmer in kmer_generator(seq.sequence, k)])
    if len(dic) < len(seq_lst):
        for s in seq_lst:
            counter += 1
            if counter % 50 == 0:
                print("sequence number: ", counter)

            all_kmers_but1 = set([kmer for seq in seq_lst if seq.name != s.name for kmer in kmer_generator(seq.sequence, k)])
            if len(all_kmers.difference(all_kmers_but1)) == 0:
                k += 1
                counter = 0
                dic= {}
                #dic[s.name] = all_kmers.difference(all_kmers_but1)
                break
            elif len(all_kmers.difference(all_kmers_but1)) > 0:
                dic[s.name] = all_kmers.difference(all_kmers_but1)
    else:
        break
        
print('Minimum length k that give unique probes for each given DNA sequences: ', k)

