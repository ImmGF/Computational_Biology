import numpy as np
from hmmlearn import hmm
from Bio import SeqIO
from collections import Counter

#funtion converts letter-like nucleotide sequences into numeric sequences read-into hmm model
def seqRecordToObservations(sequence, nucleotidePairs):
    seq = sequence.seq
    return [[o] for o in [nucleotidePairs.index(pair[0] + pair[1])
                          for pair in zip(seq[1:], seq[: -1])]]

#generate 5 different predictions
for _ in range(5):

    #define input parameters to hmm-functions from hmmlearn package:

    #- define 2 states (equivalents of 'CpG' and 'non-CpG')
    states = ["state1", "state2"]
    n_states = len(states)

    #- define two-letter nucleotide sequences which would be 'searched' by hmm models
    nucleotides = ["A", "T", "G", "C"]
    nucleotidePairs = [
        nuc1 + nuc2 for nuc1 in nucleotides for nuc2 in nucleotides]

    #-read in sequence on which hmm model would be trained on
    sequence = SeqIO.read("./cpg.fa", "fasta")

    #- define observation parameters red-into hmm model
    observations = seqRecordToObservations(sequence, nucleotidePairs)
    n_observations = len(observations)

    #- define set of probability matrixes read-into hmm model
    start_probability = np.array([0.6, 0.4])

    transition_probability = np.array([
        [0.8, 0.2],
        [0.2, 0.8]
    ])

    emission_probability = np.array([
        [0.5] * n_observations,
        [0.5] * n_observations
    ])

    #train model on cpg.fa file
    model = hmm.MultinomialHMM(n_components=n_states, verbose=True)
    model.startprob_ = start_probability
    model.transmat_ = transition_probability
    model.emissionprob_ = emission_probability

    model.fit(observations, [n_observations])

    #print transition end emission matrixes changed after training phase
    print("Transition matrix")
    print(model.transmat_)
    print("Emission probs")
    [print(f"{tup[0]}: {tup[1][0]}, {tup[1][1]}") for tup in zip(
        nucleotidePairs, zip(model.emissionprob_[0], model.emissionprob_[1]))]

    stateSequence = model.predict(observations, [n_observations])
    print(Counter(stateSequence))
    print("assume \"0\" is CpG")

    #Generate predictions on cpg_test.fa file using model trained on cpg.fa file
    import time
    with open("predictions_{0}.txt".format(time.time()), "w") as file:
        id = 0
        for seq in SeqIO.parse("./cpg_test.fa", "fasta"):
            observations = seqRecordToObservations(seq, nucleotidePairs)
            stateSequence = model.predict(observations, [len(observations)])
            counter = Counter(stateSequence)
            print('-'*40)
            print("CpG: {0}".format(counter[0]))
            print("non-CpG: {0}".format(counter[1]))
            prediction = "CpG" if counter[0] > counter[1] else "non-CpG"
            file.write(f"{seq.name} {id}: {prediction}\n")
            id += 1
