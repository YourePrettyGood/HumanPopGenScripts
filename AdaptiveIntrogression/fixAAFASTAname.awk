#!/bin/awk -f
#Renames the FASTA contig name for ancestral allele scaffolds produced by
# Ensembl, since they have a structured name, but many tools require the
# name to be an exact match to the chromosome.
#We assume here that the structure is colon-delimited and that the chromosome
# ID is found in the third element of this structure.
#Input is a FASTA of ancestral allele calls from Ensembl, e.g. as provided
# by hmmix or ArchaicSeeker2.0
#Output is a renamed FASTA that will work with e.g. bcftools +fill-from-fasta
/^>/{
   n=split($1, a, ":");
   print ">"a[3];
}
!/^>/{
   print;
}
