#!/usr/bin/env python
# coding: utf-8

# In[30]:


"""TASK 1"""

#Identify the closest related protein-coding sequences in the E. coli genome 
#for each of the input protein sequences given in the fasta file (protein_fragments.fa).

#create blast DataBase
from Bio.Blast.Applications import NcbimakeblastdbCommandline
cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file = 'data/DataBase/genes_e_coli_new.fa_', title = "blastDB")

cline()

#perform local search using blastn commandlinde
from Bio.Blast.Applications import NcbitblastnCommandline
blastx_cline = NcbitblastnCommandline(query="data/protein_fragments.fa_", db="data/DataBase/genes_e_coli_new.fa_", 
                                     evalue=0.001, outfmt=5, out="outputs/output_1/proteins_blasted.xml")

stdout, stderr = blastx_cline()

#parse the results and obtain:
# - input sequence id, best matching E. coli gene id
# - the associated e-value.
from Bio.Blast import NCBIXML

result_handle = open("outputs/output_1/proteins_blasted.xml")
blast_records = NCBIXML.parse(result_handle)

blast_dic = {}
for blast_record in list(blast_records):
    for i in range(len(blast_record.alignments)):
        protein_fragment = blast_record.alignments[i].hit_def.split()[0]
        if blast_record.query in blast_dic:
            if blast_dic[blast_record.query][1] > blast_record.descriptions[i].e:
                blast_dic[blast_record.query][0] = protein_fragment
                blast_dic[blast_record.query][1] = blast_record.descriptions[i].e
        else:
            blast_dic[blast_record.query] = []
            blast_dic[blast_record.query].append(protein_fragment)
            blast_dic[blast_record.query].append(blast_record.descriptions[i].e)

#write input sequence id, best matching E. coli gene id and the associated e-value in comma-sepparated list
def writeBestMatches(dic):
    f = open("outputs/output_1/{}.csv".format("bestMatches"), "w")
    
    with f:
        write_handle = "{},{},{}\n".format("input sequence id","best matching E. coli gene id","e-value")
        f.write(write_handle)
        for key in dic:
            write_handle = "{},{},{}\n".format(key.split()[0],dic[key][0],dic[key][1])
            f.write(write_handle)
            
writeBestMatches(blast_dic)

"""TASK 2"""

from Bio import Seq,SeqIO

#read in fasta file with promotor sequences
promoter_sequences = list(SeqIO.parse(open("data/proms_e_coli_fixed.fa_"), "fasta" ))

#group promotors' sequences to the groups (groupA and groupB) using their gene name accordingly

groupA = []
groupB = []
for blast_key in blast_dic:
    for promoter in promoter_sequences:
        if blast_key[0:6] == 'groupA':
            if blast_dic[blast_key][0] == promoter.id:
                groupA.append(promoter)
        elif blast_key[0:6] == 'groupB':
            if blast_dic[blast_key][0] == promoter.id:
                groupB.append(promoter)

#write the promoters in separate fasta files

write_handleA = open('outputs/output_2/promoters_groupA.fa', 'w')
for promoterA in groupA:
    SeqIO.write(promoterA, write_handleA, 'fasta')
    

write_handleB = open('outputs/output_2/promoters_groupB.fa', 'w')
for promoterB in groupB:
    SeqIO.write(promoterB, write_handleB, 'fasta')

from numpy import argmax
from typing import List
from tqdm import tqdm
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import Alphabet

#Calculate log-odds matrix according to consensus algorithm given motif length.
def consensus(fastaFile, motifLength, returnCount = 5):
    
    alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    seqs = [setAlphabet(seqRecord.seq, alphabet) for seqRecord in SeqIO.parse(fastaFile, "fasta")]
    firstSeq = seqs.pop(0)
    
    initialMotifs = [initialMotif for initialMotif in getSubseqsOfGivenLengths(firstSeq, motifLength)]
    
    pwms = [getPwmForAllSeqs(initialMotif,seqs) for initialMotif in initialMotifs]
    
    return(getBestScoringModels(pwms, returnCount))

def getBestScoringModels(pwms, returnCount):
    pwmsWithScores = [(pwm.mean(), pwm) for pwm in pwms]
    pwmsWithScores.sort(key=lambda x: x[0], reverse=True)
    return([pwmWithScore[1] for pwmWithScore in pwmsWithScores][:returnCount])

def setAlphabet(seq, alphabet):
    return(Seq(str(seq.upper()), alphabet = alphabet))

def getSubseqsOfGivenLengths(seq, length):
    startPoints = range(len(seq) - length)
    
    for startPoint in startPoints:
        endPoint = startPoint + length
        subSeq = seq[startPoint:endPoint]
        yield(subSeq)

def updateModel(MotifSeqs, newMotifSeq):
    MotifSeqs.append(newMotifSeq)
    return(getPwm(MotifSeqs))

#get best motif from following sequence based on log-odds of motif occurence matrix
def getBestMotifFromSeq(pwm,seq):
    scores = pwm.calculate(seq)
    bestScoreIndex = argmax(scores)
    motifLength = pwm.length
    bestSubSeq = seq[bestScoreIndex : bestScoreIndex + motifLength]
    return(bestSubSeq)

#Calculate log-odds of motif occurence matrix with heuristics following best motifs from following sequences
def getPwmForAllSeqs(InitialMotifSeq, seqs):
    motifSeqs = [InitialMotifSeq]
    pwm = getPwm(motifSeqs)
    
    for seq in seqs:
        newMotifSeq = getBestMotifFromSeq(pwm,seq)
        pwm = updateModel(motifSeqs, newMotifSeq)
        
    return(pwm)

#Calculate log-odds of a motif occurence
def getPwm(motifSeqs):
    motif = motifs.create(motifSeqs)
    counts = motif.counts.normalize(pseudocounts=1.0)
    pwm_ = counts.log_odds()
    return(pwm_)

#obtain motif profiles as pwm matrices

motifsA = consensus("outputs/output_2/promoters_groupA.fa", 15, returnCount=10)
motifsB = consensus("outputs/output_2/promoters_groupB.fa", 15, returnCount=10)

#function that writes motif profiles into .pfm format file

def writePFM(Motifs, fileName):
    
    
    alphabet = "ACGT"
    
    
    
    for i, Motif in enumerate(Motifs):
        f = open("outputs/output_2/{}_pfm/motif_{}.pfm".format(fileName,i), "w")
        with f:
            for key in alphabet:
                #print(" ".join(map(str, Motif[key])))
                f.write(" ".join(map(str, Motif[key])))
                f.write("\n")
        

# write motifs into .pfm files

writePFM(motifsA, "motifsA")
writePFM(motifsB, "motifsB")

"""TASK 3"""

#get number of positions in a given set of promoter sequences that 
#have a log-odds score higher than 0 for a given motif
def getHits(motif, group):
    
    counter = 0
    for promoter in group:
        for i in range(len(promoter.seq) - len(motif[0]) + 1):
            score = 0
            for j in range(len(motif[0])):
                score += motif[promoter.seq[i+j]][j]
            if score > 0:
                counter += 1
                
    return(counter)

#get number of combined potential positions of motif for a given group
def getPositions(group):
    all_positions = len(group) * (len(group[0]) - 15)
    return(all_positions)

from scipy.stats import binom_test
#perform a binomial test for enrichment for a given motif and get the p-value
def testBinomial(motif, group):
    groups = [groupA, groupB]
    groups.remove(group)
    alternate_group = groups[0]
    
    x = getHits(motif, group)
    n = getPositions(group)
    p = getHits(motif, alternate_group)/getPositions(alternate_group)
    
    p_value = binom_test(x, n, p, alternative='greater')
    
    return(p_value)

#write results of binomial test of enrichment in tabular format
def writeBinom(motifs, fileName, group, alternate_group):
    f = open("outputs/output_3/{}_binomial.txt".format(fileName), "w")
    
    groupName = "group{}".format(fileName[-1])
    alternate_groupName = "group{}".format("B" if fileName[-1] == "A" else "A")
    
    with f:
        write_handle = "{} \t {} \t {} \t {}\n".format("motif sequence","no hits in {}".format(groupName), "no hits in {}".format(alternate_groupName), "p-value")
        f.write(write_handle)
        for motif in motifs:
            write_handle = "{} \t {} \t {} \t {}\n".format(motif.consensus,getHits(motif, group),getHits(motif, alternate_group), testBinomial(motif, group))
            f.write(write_handle)

writeBinom(motifsA, "motifsA", groupA, groupB)

writeBinom(motifsB, "motifsB", groupB, groupA)

for motif in motifsA:
    print(motif.consensus, ": ", getHits(motif, groupA), "\t", getHits(motif, groupB),"\t", testBinomial(motif, groupA))

for motif in motifsB:
    print(motif.consensus, ": ", getHits(motif, groupB), "\t", getHits(motif, groupB),"\t", testBinomial(motif, groupB))

