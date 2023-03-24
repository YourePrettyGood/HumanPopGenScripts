#!/bin/awk -f
#This script takes a sorted version of the ArchaicSeeker2
# .seg file and outputs a BED file of the tracts with
# tract ID in the fourth column as a key-value pair.
#The output format is designed to mimic that produced
# in the sprime.nf pipeline as input for SprimeTractGeneList.awk.
#The .seg file should be sorted and passed as the first
# positional argument with a command substitution like this:
# <(tail -n+2 [.seg file] | sort -k2,2V -k3,3n -k4,4n -k1,1V | \
#   cat <(head -n1 [.seg file]) -)
#
# Input files:
#  1) ArchaicSeeker2 .seg file, sorted by region and *then* by haplotype
#
# Options:
#  pop: (required) Name of the ArchaicSeeker2.0 target population (i.e. not
#                  the outgroup)
#                  This value is only used to make a unique Tract ID
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(pop) == 0) {
      print "pop variable is missing, please set it." > "/dev/stderr";
      print "pop is the name of the ArchaicSeeker2.0 target population" > "/dev/stderr";
      print "i.e. the non-outgroup population in the ArchaicSeeker2.0 run" > "/dev/stderr";
      print "It is used as a prefix for a unique tract identifier" > "/dev/stderr";
      exit 2;
   };
}
NR==1{
   #Map column names to indices:
   for (i=1; i<=NF; i++) {
      segcols[$i]=i;
   };
}
NR>1{
   #Keep track of the tracts per chromosome and haplotype:
   hapchrom=$segcols["ID"] SUBSEP $segcols["Contig"];
   tractcounts[hapchrom]++;
   #The TractID for ArchaicSeeker2 needs to include the haplotype ID:
   tractid=pop"_"$segcols["Contig"]"_"tractcounts[hapchrom]"_"$segcols["ID"];
   #Store the start, end, tract ID, and ArchaicSeeker2 HMM state:
   tractstart=$segcols["Start(bp)"];
   tractend=$segcols["End(bp)"];
   ASstate=$segcols["BestMatchedPop"];
   #Output the tract:
   print $segcols["Contig"], tractstart, tractend, "TractID="tractid";ASSTATE="ASstate, ".", "+";
}
