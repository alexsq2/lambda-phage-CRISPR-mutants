# lambda-phage-CRISPR-mutants
Comparing CRISPR generated lambda phage mutants to the wild-type target sequence

use "protospacer mismatch finder" to match DNA sequences from a high-throughput sequencing run to a default or wild-type sequence.
you must enter your wild-type sequence, the file you wish to analyze, and the flanking sequences to extract DNA sequences from different lines in the file

use "deletion finder" to find deletions in DNA sequences from high-throughput sequencing
the script was used to find deletions in sequences generated by Pacbio sequencing
you must enter your wild-type sequence and the file you wish to analyze
the script will determine which lines contain a DNA sequence and match them to the wild-type sequence

Both scripts were run with .fastq files as inputs with DNA sequences on a single line, but other text files formats separated by lines should work
