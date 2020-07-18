#!/usr/bin/env python
# coding: utf-8

# # Assignment 1
# ## Coronavirus phylogeny


#############################
# Search for good sequences #
#############################


# access to Entrez database

from Bio import SeqIO
from Bio import Entrez
Entrez.email="if290707@mimuw.edu.pl"
Entrez.tool="biopython"


# 0. Create semi-automatized dictionary of expressions 'Search_terms_dictionary' and 
#    function 'entrez_term_creator' whose aim is to narrow down query results for sequences
#    of searched virus.

#    Organism:
#    *2019-nCov virus (current coronavirus)
#    *SARS – another coronavirus back from 2002
#    *bat coronavirus
#    *MERS – middle eastern respiratory syndrome coronavirus
#    *influenza A -  flu virus, to see how far are coronaviruses from the flu
#    *hepatitis A virus – less related virus


# 'Search_terms_dictionary' keys:
#  - [virus name] - name of the virus searched in the DB
#  - [protein term] - terms for proteins searched along with the virus name
#  - [host] - name of the host of the searched virus
#  - [exclusion comments] - terms that exclude given query to be the output (e.g. 
#    if in the searched query there are 'Severe acute respiratory syndrome-related coronavirus' and 'Wuhan'
#    then query is not qualified as viable query result)

Search_terms_dictionary = {
    "virus name": ["Severe acute respiratory syndrome coronavirus 2", 
                   "Severe acute respiratory syndrome-related coronavirus",
                   "Bat SARS-like coronavirus",
                   "Middle East respiratory syndrome-related coronavirus",
                   "Influenza A virus",
                   "Hepatovirus A"
                  ],
    "protein term": [["spike glycoprotein", "spike protein"],
                     ["spike glycoprotein", "spike protein"],
                     ["spike glycoprotein", "spike protein"],
                     ["spike glycoprotein", "spike protein"],
                     ["glycoprotein", "hemagglutinin", "neuraminidase"],
                     ["capsid protein"]
                    ],
    "host": [["Homo sapiens"],
             ["Homo sapiens"],
             ["Rhinolophus sinicus"],
             ["Homo sapiens"],
             [],
             ["Homo sapiens"]
            ],
    "exclusion comments": [[],
                          ["Wuhan","Severe acute respiratory syndrome coronavirus 2"],
                          [],
                          [],
                          ["cds"],
                          ["cds"]]
}




# Function 'entrez_term_creator' generates query terms based on the 'Search_terms_dictionary'.
# The function generates separate terms for spike protein query and virus genome query.

def entrez_term_creator():
    protein_terms = []
    genome_terms = []
    
    for i in range(6):
        term = ""
        term += Search_terms_dictionary["virus name"][i] + " [ORGN] "
        if Search_terms_dictionary["protein term"][i]:
            term +=  " AND " + "(" + " OR ".join([ p  + " [PROT] " for p in Search_terms_dictionary["protein term"][i]]) + ")"
        if Search_terms_dictionary["host"][i]:
            term += " AND " + "host " + Search_terms_dictionary["host"][i][0] + " [TEXT] "
        if Search_terms_dictionary["exclusion comments"][i]:
            term += " NOT " + "("+ " AND ".join([e + " [ALL]" for e in Search_terms_dictionary["exclusion comments"][i]]) + ")"
        protein_terms.append(term)
    
    for j in range(6):
        term = ""
        term += Search_terms_dictionary["virus name"][j] + " [ORGN] "
        if Search_terms_dictionary["host"][j]:
            term += " AND " + "host " + Search_terms_dictionary["host"][j][0] + " [TEXT] "
        term += " AND whole genome[ALL]"
        if Search_terms_dictionary["exclusion comments"][j]:
            term += " NOT " + " OR ".join([e + " [ALL]" for e in Search_terms_dictionary["exclusion comments"][j]])
        
        
        genome_terms.append(term)
    
    return(genome_terms,protein_terms)



gen_terms, prot_terms = entrez_term_creator()



# This piece of short loop is controlled manually. This step requires manual setup of 'virus_idx' (virus index
# corresponds to the order of virus from Search_terms_dictionary["virus name"]). As the result loop prints out
# a number of queries (i.e. 15), which the programmer would manually choose as the best fit.
# The loop generates Entrez id of the query, which then should be fed into proteins_to_align and genomes_to_align
# lists which keep name of the virus and its' chosen Entrez id (which will be directly used to download 
# sequence of itnerest).

# Loop for spike proteins:

virus_idx = 0
print("#####################")
print(Search_terms_dictionary["virus name"][virus_idx])
handle=Entrez.esearch(db="protein",term= prot_terms[virus_idx], retmax=15)
spike_rec=Entrez.read(handle)
for i in spike_rec["IdList"]:
    print("==================")
    print("protein id:", i)
    spike_handle=Entrez.efetch(db="protein",id=i,rettype="fasta")
    spike_virus=SeqIO.read(spike_handle,"fasta")
    print(spike_virus.description)
    print("długość: ", len(spike_virus.seq))
    print("==================")


print("#####################")


# Loop for viral genomes:

virus_idx = 4
print("#####################")
print(Search_terms_dictionary["virus name"][virus_idx])
handle=Entrez.esearch(db="nuccore",term= "Hepatovirus A [ORGN]  AND genome[ALL] NOT (segment[ALL] OR protein[ALL] OR partial[ALL] OR gene[ALL])", retmax=30)
spike_rec=Entrez.read(handle)
for i in spike_rec["IdList"]:
    print("==================")
    print("genome id:", i)
    spike_handle=Entrez.efetch(db="nucleotide",id=i,rettype="fasta")
    spike_virus=SeqIO.read(spike_handle,"fasta")
    print(spike_virus.description)
    print("długość: ", len(spike_virus.seq))
    print("==================")


print("#####################")


# *As search queries using expressions regarding Influenza A and Hepatovirus A proteins both of those virus
#  do not have spike proteins but posses coat proteins which were used. 
#  For Influenza A two proteins were used:
#  - neuraminidase
#  - hemagglutinin
#  For Hepatovirus A:
#  - short sequence capsid
#  - long sequence capsid
#  - long sequence coat protein not identifies as capsid protein

proteins_to_align = [("2019-nCov virus", 1807246654),
                     ("SARS 2002", 30795145),
                     ("Bat SARS", 1369125431),
                     ("MERS", 1783577429),
                     ("Flu A neuraminidase", 1802172444),
                     ("Flu A hemagglutinin", 1802172426),
                     ("Hepatovirus A big protein", 1805598189),
                     ("Hepatovirus A big capsid ", 14029755),
                     ("Hepatovirus A small capsid", 6449363)
                    ]
genomes_to_align = [("2019-nCov virus", 1804153872),
                    ("SARS 2002", 32493129),
                    ("Bat SARS", 1802633852),
                    ("MERS", 1707416314),
                    ("Flu A virus", 92582873),
                    ("Hepatovirus A", 1805598188),
                   ]


from Bio import SeqIO, Seq, Phylo, AlignIO

#1a. Downloading the spike protein sequences from Entrez using Bio.Entrez

proteins_to_fasta = []
genomes_to_fasta = []

for prot in proteins_to_align:
    
    prot_handle=Entrez.efetch(db="protein",id=str(prot[1]),rettype="fasta")
    
    virus_protein=SeqIO.read(prot_handle,"fasta")
    virus_protein.id = prot[0].replace(" ", "_")
    
    proteins_to_fasta.append(virus_protein)
    
#Write chosen viral protein sequences to fasta file
with open("virus_spike_proteins.fasta", "w") as text_handle:
    SeqIO.write(proteins_to_fasta, text_handle, "fasta")

#1b. Downloading the viral genome sequences from Entrez using Bio.Entrez
    
for genome in genomes_to_align:
    
    genome_handle=Entrez.efetch(db="nucleotide",id=str(genome[1]),rettype="fasta")
    
    virus_genome=SeqIO.read(genome_handle,"fasta")
    virus_genome.id = genome[0].replace(" ", "_")
    
    genomes_to_fasta.append(virus_genome)
    
#Write chosen viral genome sequences to fasta file
with open("virus_genomes.fasta", "w") as text_handle:
    SeqIO.write(genomes_to_fasta, text_handle, "fasta")


from Bio.Align.Applications import ClustalwCommandline

#2. Creating multiple alignments of both DNA (whole genome) and 
#   protein (spike protein) sequences

# set clustal handle
handler_clustal = '/usr/bin/clustalw'

#run multiple alignment algorithm (ClustalW) on DNA and protein sequences on given virus
clustal_spike_alignment = ClustalwCommandline(handler_clustal, infile='virus_spike_proteins.fasta')
clustal_genome_alignment = ClustalwCommandline(handler_clustal, infile='virus_genomes.fasta')

out_log, err_log = clustal_spike_alignment()
out_log, err_log = clustal_genome_alignment()


from Bio.Phylo.TreeConstruction import *
#3. Creating phylogenetic trees for both protein and 
#   genome alignments using the UPGMA and neighbor joining algorithms

spike_aln = AlignIO.read('virus_spike_proteins.aln', 'clustal')
genome_aln = AlignIO.read('virus_genomes.aln', 'clustal')

constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')

spike_dm = calculator.get_distance(spike_aln)
genome_dm = calculator.get_distance(genome_aln)

constructor.upgma
#UPGAMA
upgma_spiketree = constructor.upgma(spike_dm)
upgma_genometree = constructor.upgma(genome_dm)

#NJ
nj_spiketree = constructor.nj(spike_dm)
nj_genometree = constructor.nj(genome_dm)


from matplotlib import pyplot as plt

#4. Visualizing the trees
#UPGAMA
#Phylo.draw_ascii(upgma_spiketree)
#Phylo.draw_ascii(upgma_genometree)
plt.rcParams["figure.figsize"] = [15, 8]
Phylo.draw(upgma_spiketree)
Phylo.draw(upgma_genometree)


#NJ
#Phylo.draw_ascii(nj_spiketree)
#Phylo.draw_ascii(nj_genometree)
plt.rcParams["figure.figsize"] = [15, 8]
Phylo.draw(nj_spiketree)
Phylo.draw(nj_genometree)

