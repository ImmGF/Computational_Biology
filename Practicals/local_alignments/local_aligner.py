#!/usr/bin/env python
# coding: utf-8


from Bio import Seq,SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo

import numpy as np
import pandas as pd


def DNA_translator(sequence):
    
    protein_sequences = []
    s = Seq.Seq(sequence)
    
    for i in range(0,3):
        protein_sequences.append(s[i::].translate())
    
    return(protein_sequences)


def amino_substitution_cost(amino1, amino2):
    
    substitution_matrix = MatrixInfo.blosum60
    
    if (amino1,amino2) in substitution_matrix.keys():
        cost = substitution_matrix[(amino1,amino2)]
    elif (amino2,amino1) in substitution_matrix.keys():
        cost = substitution_matrix[(amino2,amino1)]
    
    return(cost)


def create_matrix(sequence1, sequence2, indel_cost):
    
    cols = ['#'] + [a for a in sequence2]
    rows = ['#'] + [b for b in sequence1]
    
    alignment_matrix = pd.DataFrame(np.zeros((len(sequence1) + 1, len(sequence2) + 1), 'int'))
    
    transition_matrix = pd.DataFrame(np.zeros((len(sequence1) + 1, len(sequence2) + 1), 'str'))
    
    for s1 in range(0,len(sequence1)):
        for s2 in range(0,len(sequence2)):
            
            if sequence1[s1] != '*' and sequence2[s2] != '*':
                
                alignment_matrix.iloc[s1+1,s2+1] = max(0,
                                             alignment_matrix.iloc[s1+1,s2] - indel_cost,
                                             alignment_matrix.iloc[s1,s2+1] - indel_cost,
                                             alignment_matrix.iloc[s1,s2] + amino_substitution_cost(sequence1[s1],sequence2[s2]))
                if alignment_matrix.iloc[s1+1,s2+1] != 0:
                    if alignment_matrix.iloc[s1+1,s2+1] == alignment_matrix.iloc[s1+1,s2] - indel_cost:
                        transition_matrix.iloc[s1+1,s2+1] = 'lft'
                    elif alignment_matrix.iloc[s1+1,s2+1] == alignment_matrix.iloc[s1,s2+1] - indel_cost:
                        transition_matrix.iloc[s1+1,s2+1] = 'up'
                    elif alignment_matrix.iloc[s1+1,s2+1] == alignment_matrix.iloc[s1,s2] + amino_substitution_cost(sequence1[s1],sequence2[s2]):
                        transition_matrix.iloc[s1+1,s2+1] = 'diag'
    
    for s1 in range(0,len(sequence1)):
        for s2 in range(0,len(sequence2)):
    
            if sequence1[s1] == '*' :
                
                alignment_matrix.iloc[s1+1,:] = -float('inf')
                transition_matrix.iloc[s1+1,:] = 'inf'
                
            elif sequence2[s2] == '*':
                
                alignment_matrix.iloc[:,s2+1] = -float('inf')
                transition_matrix.iloc[:,s2+1] = 'inf'
    
    return(alignment_matrix, transition_matrix)


def align_sequences(alg_matrix, trn_matrix, sequence1, sequence2):
    
    local_alignment_score = alg_matrix.max(1).max()
    
    col = alg_matrix.max(0).idxmax() - 1
    row = alg_matrix.max(1).idxmax() - 1
    
    sequence1_end = row
    sequence2_end = col
    
    direction = trn_matrix.loc[row+1,col+1]
    
    aligned_sequence1 = ''
    aligned_sequence2 = ''
    
    while direction != '' and direction != 'inf':
        
        if direction == 'diag':
            aligned_sequence1 += sequence1[row]
            aligned_sequence2 += sequence2[col]
            row -= 1
            col -= 1
            
        elif direction == 'lft':
            aligned_sequence1 += '-'
            aligned_sequence2 += sequence2[col]
            row -= 0
            col -= 1
        elif direction == 'up':
            aligned_sequence1 += sequence1[row]
            aligned_sequence2 += '-'
            row -= 1
            col -= 0
        direction = trn_matrix.loc[row + 1,col + 1]
        
    sequence1_start = row + 1
    sequence2_start = col + 1
    
    return(local_alignment_score, aligned_sequence1[::-1],aligned_sequence2[::-1], [sequence1_start,sequence1_end], [sequence2_start,sequence2_end])


def align_codons(dna_seq1, dna_seq2, codon_deletion_cost):
    
    protein_seqs1 = DNA_translator(dna_seq1)
    protein_seqs2 = DNA_translator(dna_seq2)
    
    score_matrix = pd.DataFrame(np.zeros((len(protein_seqs1), len(protein_seqs2)), 'float'))
    coordinates1_matrix = {}
    coordinates2_matrix = {}
    aligned_protein1_matrix = {}
    aligned_protein2_matrix = {}
    
    for prot_seq1_no in range(len(protein_seqs1)):
        
        coordinates1_matrix[prot_seq1_no] = []
        coordinates2_matrix[prot_seq1_no] = []
        aligned_protein1_matrix[prot_seq1_no] = []
        aligned_protein2_matrix[prot_seq1_no] = []
        
        for prot_seq2_no in range(len(protein_seqs2)):
            
            alg, trn = create_matrix(protein_seqs1[prot_seq1_no], protein_seqs2[prot_seq2_no], codon_deletion_cost)
            score, aligned_protein1, aligned_protein2, coordinates1, coordinates2 = align_sequences(alg, trn, protein_seqs1[prot_seq1_no], protein_seqs2[prot_seq2_no])
            
            score_matrix.iloc[prot_seq1_no,prot_seq2_no] = score
            
            aligned_protein1_matrix[prot_seq1_no].append(aligned_protein1)
            aligned_protein2_matrix[prot_seq1_no].append(aligned_protein2)
            
            coordinates1_matrix[prot_seq1_no].append(coordinates1)
            coordinates2_matrix[prot_seq1_no].append(coordinates2)
            
            """
            print("XXXXX")
            print("DNA1_frame: ", prot_seq1_no, "DNA2_frame: ",prot_seq2_no)
            print("Score: ", score)
            print(aligned_protein1)
            print(aligned_protein2)
            print("XXXXX")
            """
    
    key = score_matrix.max(1).idxmax()
    lst_element = score_matrix.max(0).idxmax()
    
    final_score = score_matrix.iloc[key, lst_element]
    cr1 = coordinates1_matrix[key][lst_element]
    cr2 = coordinates2_matrix[key][lst_element]
    
    
    aligned_and_translated_dna_seq1 = dna_seq1[key::][3 * cr1[0] : 3 * (cr1[1] + 1)]
    aligned_and_translated_dna_seq2 = dna_seq2[lst_element::][3 * cr2[0] : 3 * (cr2[1] + 1)]
    
    dna_aligned_sequence1 = ""
    dna_aligned_sequence2 = ""
    
    counter_i = -1
    
    for i in range(len(aligned_protein1_matrix[key][lst_element])):
        if aligned_protein1_matrix[key][lst_element][i] != '-':
            counter_i += 1
            dna_aligned_sequence1 += "(" + aligned_and_translated_dna_seq1[3*counter_i:3*counter_i + 3] + ")"
            
        elif aligned_protein1_matrix[key][lst_element][i] == '-':
            dna_aligned_sequence1 += "(" + 3*'-' + ")"
    
    counter_j = -1
    
    for j in range(len(aligned_protein2_matrix[key][lst_element])):
        
        if aligned_protein2_matrix[key][lst_element][j] != '-':
            counter_j += 1
            dna_aligned_sequence2 += "(" + aligned_and_translated_dna_seq2[3*counter_j:3*counter_j + 3] + ")"
            
        elif aligned_protein2_matrix[key][lst_element][j] == '-':
            dna_aligned_sequence2 += "(" + 3*'-' + ")"
    
    """
    print("############")
    
    print(aligned_protein1_matrix[key][lst_element])
    print(aligned_protein2_matrix[key][lst_element])
    
    print(protein_seqs1[key][cr1[0]:cr1[1]+1])
    print(protein_seqs2[lst_element][cr2[0]:cr2[1]+1])
    
    print()
    """
    
    return(final_score,dna_aligned_sequence1,dna_aligned_sequence2, aligned_protein1_matrix[key][lst_element], aligned_protein2_matrix[key][lst_element])


DNA1 = 'AGGTGGTAATAAACTAGTGGGTAATTAAAGACTTTTAGTAGGATCTAATGCTAAGGCTGGGATAATCTGCTAAGGCTGGGATAATCTGCTAAGGCTGGTTTTAAGACATTCCCAGCCGTGCACAGATTGCTTAAAAGAAG'
DNA2 = 'ATGTCCGGTAATGGTAAAGTAAGTGGTAAATAAGCTGGTTTAACAGCTGCTAAAGCTTTAACTCAATCTAGATCTGCTAAGGCTTAAGGTTTGACATTCCCAGTCGGTATAAGAGTGCACAGATTGCTAAGAATAAGAGGTAAC'


score, aligned_dna_sequence_no1, aligned_dna_sequence_no2, aligned_protein_no1, aligned_protein_no2 = align_codons(DNA1, DNA2, 1)


print("Score: ", score)
print()
print("Aligned protein sequences:")
print(aligned_protein_no1)
print(aligned_protein_no2)
print()
print("Aligned DNA sequences by codons:")
print(aligned_dna_sequence_no1)
print(aligned_dna_sequence_no2)

