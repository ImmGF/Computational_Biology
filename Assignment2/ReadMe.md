Try to combine a few of the problems we have encountered: Searching for related sequences, identifying regulatory regions and motifs and enrichment analysis. All of This will be performed on a provided dataset and the expected result will be a program that solves the problem and a report with your findings.
 
1. Identify the closest related protein-coding sequences in the E. coli genome for each of the input protein sequences given 
in the fasta file (protein_fragments.fa). This should be done by downloading the E. coli protein-coding sequences 
(from this file – genes_e_coli_new.fa), creating a BLAST database from them and performing a local BLAST search against 
this database. Importantly, we have protein fragments that need to be matched against DNA sequences, so we have to do something
to make this compatible. (Here’s a tutorial on running blast locally.) The result of this step should be a comma-separated list 
containing: input sequence id, best matching E. coli gene id and the associated e-value.

2. For each of the identified E. coli genes, you should find the associated promoter DNA sequence in the given file (proms_e_coli_new.fa, 
proms_e_coli_fixed.fa). Please note that the input file contained sequences from groups A and B,  we speculate (based on the 
empirical evidence), that these groups of genes should have different regulatory mechanisms. We need to identify 10 sequence 
motifs present in the promoters associated with group A and 10 motifs associated with group B, independently of each other. 
All motifs should be of length 15 . The result of this step should be two sets of 10 different motif position-specific matrices 
in a .pfm format. You can do this with a consensus-like method you implement yourself, or you can use the MEME-suite and
process its output.

3. Given the two sets of motifs, we would like to select only the motifs that are specific to group A or group B. For this, we need to compute the occurrences of the motifs from .pfm in promoters from both groups. For this purpose, we should scan the promoter sequences of genes from groups A and B, and obtain a number of positions in these promoter sequences that have a log-odds score higher than 0. Now, given the combined lengths of the promoters in each group and the total number of “hits” in each group, we can use the binomial test to identify the motifs that are significantly enriched in one of the groups against in the other. The result of this step should be a list of all 20 motifs with associated, number of hits in the promoters from group A, and promoters from group B and the associated p-values for enrichment in group A and group B respectively.

All of these results and your implementation, should be described in a short report.
