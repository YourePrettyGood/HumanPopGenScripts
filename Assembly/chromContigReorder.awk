#!/bin/awk -f
#This script takes FASTA indices (.fai) of the query contigs and reference
# chromosomes as well as the sorted output of gpBestHitList.awk and outputs
# new chromosome and contig orders for mummerplot on STDERR and STDOUT,
# respectively. The files produced on STDERR and STDOUT can be used with
# the -R and -Q options of mummerplot, respectively.
#Inputs:
# [query contigs/scaffolds FAI]
# [reference chromosomes FAI]
# gpBestHitList.awk [-fat .gp] [BEDbestHit.awk output] | sort -k2,2V -k3,3n
#So a full run might look like:
#nucmer -l100 -c1000 -t16 -p [query]_[ref] [query FASTA] [ref FASTA]
#mummerplot --png --large --fat --filter -p [query]_[ref]_l100c1000_filtered [query]_[ref].delta
#show-coords -cdHloqT [query]_[ref].delta | \
# coordToBED.awk | \
# sort -k4,4V -k1,1V -k2,2n -k3,3n | \
# BEDbestHit.awk | \
# gpBestHitList.awk [query]_[ref]_l100c1000_filtered.gp - | \
# sort -k2,2V -k3,3n | \
# chromContigReorder.awk [query FASTA].fai [ref FASTA].fai - 2> [query]_vs_[ref]_chrom_order.tsv > [query]_vs_[ref]_ctg_order.tsv
#mummerplot --png --large -Q [query]_vs_[ref]_ctg_order.tsv -R [query]_vs_[ref]_chrom_order.tsv --filter -p [query]_[ref]_l100c1000_filtered_reordered [query]_[ref].delta
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#.fai of the query contigs/scaffolds:
#We store a map of contigs/scaffolds to their lengths.
filenum==1{
   querylen[$1]=$2;
}
#.fai of the reference chromosomes:
#We store a map of chromosomes to their lengths.
filenum==2{
   reflen[$1]=$2;
}
#
filenum==3{
   #Set orientation of query contig/scaffold:
   print $1, querylen[$1], $5;
   #Set orientation of ref chromosome to + and output in the order
   # of the contig mappings:
   if (!($2 in refs)) {
      print $2, reflen[$2], "+" > "/dev/stderr";
      refs[$2]=1;
   };
}
