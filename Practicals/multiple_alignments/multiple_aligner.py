#!/usr/bin/env python
# coding: utf-8

# # MULTIPLE ALIGNER


from Bio import Seq,SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo

histones = list(SeqIO.parse(open("histones.fa"),"fasta"))
bzips = list(SeqIO.parse(open("bzips.fa"),"fasta"))


#Function creates a list of sequences and names of the leaves which will be during alignment used to create 
#tree-like alignment

def create_leaves(sequences):
    
    leaves = [(sequences[i], "s{}".format(i)) for i in range(len(sequences))]
    
    return(leaves)


#Function compares original sequence and aligned sequence (with additional dashes '-'), so that this 
#'dash-profile' can be propagated onto other sequences under the same node in the tree-like structure

def get_gap_position(original_sequence, aligned_sequence):
    
    gaps_first = []
    i_original = 0
    i_aligned = 0
    #get dashes from the beggning of the sequence
    while original_sequence[i_original] == '-':
        i_original += 1
    
    while aligned_sequence[i_aligned] == '-':
        i_aligned += 1
    
    i = i_aligned - i_original
    gaps_first = i * [0]
    
    original_sequence_reversed = original_sequence[::-1]
    aligned_sequence_reversed = aligned_sequence[::-1]
    
    gaps_last = []
    j_original = 0
    j_aligned = 0
    #get dashes from the end of the sequence
    while original_sequence_reversed[j_original] == '-':
        j_original += 1
    
    while aligned_sequence_reversed[j_aligned] == '-':
        j_aligned += 1
        
    j = j_aligned - j_original
    gaps_last = j * [float('inf')]
    
    if j == 0:
        j = 1
    
    #get dashes from the within the sequence
    aligned_sequence = aligned_sequence[i:-j]
    
    orig_counter = -1
    align_counter = -1
    gaps_position = []
    
    for k in range(len(aligned_sequence)):

        orig_counter += 1
        align_counter += 1
        
        if original_sequence[orig_counter] != aligned_sequence[align_counter]:
            gaps_position.append(orig_counter)
            orig_counter -= 1
    
    return(gaps_first, gaps_position, gaps_last)

#Funtion applies profile of dashes produced by the get_gap_position() function onto a give sequence

def update_gaps(sequence, gaps):
    
    for gap in gaps[1][::-1]:
        sequence = sequence[0:gap] + '-' + sequence[gap::]
        
    for gap in gaps[0]:
        sequence = '-' + sequence
        
    for gap in gaps[2]:
        sequence = sequence + '-'
        
    return(sequence)


#Function checks if alignment reached the root of the tree-like structure whose leaves-sequences are being
#alinged. If the root is reached then the function returns 'True'

def check_if_root(sequences):
    
    nodes = []
    
    for sequence in sequences:
        nodes.append(sequence[1])
        
    if len(set(nodes)) == 1:
        return(True)
    else:
        return(False)


#The funtion generates multiple alignment profile. It uses UPGMA-like greedy approach.
#1. Seeks for the best pair of the alignment within the given sequences.
#2. The pair is grouped under common node.
#3. Seeks for another best-fit pair:
#   - if best fit is paired with already aligned pair it clusters the sequence to the already existing node,
#     creating branch
#   - otherwise it creates independent pair under one node
#4. Eventually all pairs, nodes, independent sequences would be grouped under one node-root, thus creating
#   multiple alignment.

def multiple_alignment(sequences):
    
    nodes = create_leaves(sequences)
    
    while not check_if_root(nodes):
    #for x in range(4):
        
        nodes_dummy = [leaf for leaf in nodes]
        
        #Seek the best pair from within the given sequences
        maximum = -float('inf')
        aligned_pair = []
        nodes_to_replace = []
        
        for n1 in range(len(nodes)):
            nodes_without_n1 = nodes[0:n1] + nodes[n1+1:]
            for n2 in range(len(nodes_without_n1)):
                
                leaf1 = nodes[n1]
                leaf2 = nodes_without_n1[n2]
                
                
                if leaf1[1] != leaf2[1]:
                    
                    alignment = pairwise2.align.globalxx(leaf1[0], leaf2[0])
                    
                    if alignment[0][3] > maximum:
                        maximum = alignment[0][3]
                        leaf_to_write1 = (alignment[0][0],leaf1[1] + '_' + leaf2[1])
                        leaf_to_write2 = (alignment[0][1],leaf1[1] + '_' + leaf2[1])
                        idx1 = nodes.index(leaf1)
                        idx2 = nodes.index(leaf2)
        
        #Get two new alignments, oroginal sequences that were aligned and dash-profile which is to be updated
        aligned_pair = [leaf_to_write1, leaf_to_write2]
        nodes_to_replace = [idx1, idx2]
        
        gaps_of_pair = [get_gap_position(nodes[idx1][0], aligned_pair[0][0]), 
                        get_gap_position(nodes[idx2][0], aligned_pair[1][0])]
        
        #update all nodes according to the new profile from the newly created alignment
        for i in range(2):
            
            aligned_pair_element = aligned_pair[i]
            idx = nodes_to_replace[i]
            g = gaps_of_pair[i]
            
            for node_no in range(len(nodes)):
                
                if nodes[node_no][1] == nodes[idx][1]:
                    
                    new_seq = update_gaps(nodes[node_no][0], g)
                    
                    nodes_dummy[node_no] = (new_seq, aligned_pair_element[1])
                    
        nodes = nodes_dummy
    
    return(nodes)
            


sequence_1 = "AAAAAAAAAAAAAATTTTAAAATTTTCCG"
sequence_2 = "GGAATTAAAAAAAAAATTTTAAAATTTTCCG"
sequence_3 = "AAAAAAAAAAAAAAATTTTAAAATTTTCCGTTTTG"
sequence_4 = "AAAAAAAAAAAAATTAAAATTTTCCG"
sequence_5 = "AAAAAAAAAAAATTTTAAAATTTTCCGC"
sequence_6 = "CGCGCGTAAAAAAAAAAAATTTTAAAATTTTCCGC"
sequence_7 = "TTTTAAAATTTCCGC"
sequence_8 = "AATTTTAAATTTTCCGC"

unaligned_sequences = [sequence_3, sequence_4, sequence_5, sequence_1, sequence_2, sequence_6, sequence_7, sequence_8]



mlt_algn = multiple_alignment(unaligned_sequences)



for un in unaligned_sequences:
    print(un)



for ma in mlt_algn:
    print(ma[0])



histones_lst = [str(hist.seq) for hist in histones]



hist_algn = multiple_alignment(histones_lst)


for h in hist_algn:
    print(h[0][:105])

