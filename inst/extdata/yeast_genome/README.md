# File content

## yeast_cds_2001.tab

Coordinates of coding sequences in the genome of the budding yeast
(*Saccharomyces cerevisiae*), from the 2001 annotations. These
annotations were essentially done by applying an automatic rule: any
ORF of length >= 300 was considered as a putative coding gene.

## yeast_cds_2014.tab

Coordinates of coding sequences in the genome of the budding yeast
(*Saccharomyces cerevisiae*), from the 2014 annotations.

These annotations were cleaned with the cumulated results of 18 years
of additional research on the yeast genome, including transcriptome
analysis and comparative genomics. In particular, ~500 pseudo-ORFs
have been removed from the annotations.

## yeast_3nt_freq_whole_genome.tab

Trinucleotide frequencies counted in the full genome of the yeast
Sacccharomyces cerevisiae, on both strands separately, all occurrences
counted (renewing and overlapping).

## rand_genome_cds.tab

ORF predictions from a random sequence generated according to a Markov
model of order 2, whose parameters were estimated from 3-mer
frequencies in the full genome of Saccharomyces cerevisiae.

Note: I temporarily copy twice the results with 8 chromosomes, because the NCBI ORF finder does not accept sequences >1Mb, and I thus have to run it separately for each chromosome. This means that the data is a duplication of an experiment with half the genome size, but since the goal is to estimate ORF length distributions, it should be a reasonable approximation.


