#!/usr/bin/env python
# coding: utf-8


from numpy import argmax
from typing import List
from tqdm import tqdm
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import Alphabet


#Calculate log-odds of a motif occurence
def getPwm(motifSeqs):
    motif = motifs.create(motifSeqs)
    counts = motif.counts.normalize()
    pwm = counts.log_odds()
    return(pwm)


#Calculate log-odds of motif occurence matrix with heuristics following best motifs from following sequences
def getPwmForAllSeqs(InitialMotifSeq, seqs):
    motifSeqs = [InitialMotifSeq]
    pwm = getPwm(motifSeqs)
    
    for seq in seqs:
        newMotifSeq = getBestMotifFromSeq(pwm,seq)
        pwm = updateModel(motifSeqs, newMotifSeq)
        
    return(pwm)


#get best motif from following sequence based on log-odds of motif occurence matrix
def getBestMotifFromSeq(pwm,seq):
    scores = pwm.calculate(seq)
    bestScoreIndex = argmax(scores)
    motifLength = pwm.length
    bestSubSeq = seq[bestScoreIndex : bestScoreIndex + motifLength]
    return(bestSubSeq)


def updateModel(MotifSeqs, newMotifSeq):
    MotifSeqs.append(newMotifSeq)
    return(getPwm(MotifSeqs))


def getSubseqsOfGivenLengths(seq, length):
    startPoints = range(len(seq) - length)
    
    for startPoint in startPoints:
        endPoint = startPoint + length
        subSeq = seq[startPoint:endPoint]
        yield(subSeq)


def setAlphabet(seq, alphabet):
    return(Seq(str(seq.upper()), alphabet = alphabet))


def getBestScoringModels(pwms, returnCount):
    pwmsWithScores = [(pwm.mean(), pwm) for pwm in pwms]
    pwmsWithScores.sort(key=lambda x: x[0], reverse=True)
    return([pwmWithScore[1] for pwmWithScore in pwmsWithScores][:returnCount])


#Calculate log-odds matrix according to consensus algorithm given motif length.
def consensus(fastaFile, motifLength, returnCount = 5):
    
    alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    seqs = [setAlphabet(seqRecord.seq, alphabet) for seqRecord in SeqIO.parse(fastaFile, "fasta")]
    firstSeq = seqs.pop(0)
    
    initialMotifs = [initialMotif for initialMotif in getSubseqsOfGivenLengths(firstSeq, motifLength)]
    
    pwms = [getPwmForAllSeqs(initialMotif,seqs) for initialMotif in initialMotifs]
    
    return(getBestScoringModels(pwms, returnCount))



motif_length3 = consensus("ecoli_proms.fa_", 3)



motif_length4 = consensus("ecoli_proms.fa_", 4)



motif_length5 = consensus("ecoli_proms.fa_", 5)



motif_length6 = consensus("ecoli_proms.fa_", 6)



motif_length3



motif_length4


motif_length5


motif_length6

