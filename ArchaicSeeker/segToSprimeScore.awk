#!/bin/awk -f
#This script takes a sorted version of the ArchaicSeeker2
# .seg file as its first input and a pared down and subsampled
# version of the target population's VCF to produce a
# variant of the Sprime .score format.
#Keep in mind that Sprime's tracts are identified at the
# population level, but the sites output by this script are
# quantized at the haplotype level, so there are more of
# them and their records are interdigitated in the output.
#The .seg file should be sorted and passed as the first
# positional argument with a command substitution like this:
# <(tail -n+2 [.seg file] | sort -k2,2V -k3,3n -k4,4n -k1,1V | \
#   cat <(head -n1 [.seg file]) -)
#
# Required arguments:
#  chrom:    Which chromosome to focus on (since Sprime .score files
#             are per-chromosome, but .seg is for all autosomes)
#
# Input files:
#  1) ArchaicSeeker2 .seg file, sorted by region and *then* by haplotype
#  2) Output of: bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'
#      for the samples in the target population
BEGIN{
   FS="\t";
   OFS=FS;
   #Check for the required arguments:
   if (length(chrom) == 0) {
      print "Missing argument chrom. Quitting." > "/dev/stderr";
      exit 3;
   } else {
      split(chrom, chromarr, ",");
      for (i in chromarr) {
         chroms[chromarr[i]]=i;
      };
   };
   #Keep track of the file we're on:
   filenum=0;
   #Set array traversal order:
   PROCINFO["sorted_in"]="@val_num_asc";
}
FNR==1{
   #Keep track of the file we're on:
   filenum++;
}
#First file is the ArchaicSeeker2 .seg file:
filenum==1&&FNR==1{
   #Map column names to indices:
   for (i=1; i<=NF; i++) {
      segcols[$i]=i;
   };
}
filenum==1&&FNR>1{
   #Only include tracts from the chromosome(s) of interest:
   if ($segcols["Contig"] in chroms) {
      #Keep track of the tracts per chromosome and haplotype:
      hapchrom=$segcols["ID"] SUBSEP $segcols["Contig"];
      tractindex[hapchrom]=1;
      tractcounts[hapchrom]++;
      #The TractID for ArchaicSeeker2 needs to include the haplotype ID:
      tractid=$segcols["Contig"]"_"tractcounts[hapchrom]"_"$segcols["ID"];
      #We use a slightly simpler key that's still unique to each tract:
      tractkey=hapchrom SUBSEP tractcounts[hapchrom];
      #Store the start, end, tract ID, and ArchaicSeeker2 HMM state:
      tractstarts[tractkey]=$segcols["Start(bp)"];
      tractends[tractkey]=$segcols["End(bp)"];
      tractids[tractkey]=tractid;
      ASstates[tractkey]=$segcols["BestMatchedPop"];
      if (length(debug) > 0 && debug >= 2) {
         print hapchrom, tractindex[hapchrom], tractcounts[hapchrom] > "/dev/stderr";
         print tractkey, tractstarts[tractkey], tractends[tractkey], tractids[tractkey], ASstates[tractkey] > "/dev/stderr";
      };
   };
}
#Second file is the bcftools query output:
filenum==2&&FNR==1{
   #Map column names to indices:
   for (i=1; i<=NF; i++) {
      #The entire header line gets prefixed with "# ":
      sub("#[ ]?", "", $i);
      #Each column name is prefixed with it's 1-based index in brackets:
      gsub("[[0-9]+]", "", $i);
      if ($i ~ /:GT$/) {
         #Each FORMAT field is suffixed with ":[tag]":
         gsub(":GT", "", $i);
         #Keep track of the samples for iteration:
         gtcols[$i]=i;
         if (length(debug) > 0 && debug >= 1) {
            print i, $i > "/dev/stderr";
         };
      };
      querycols[$i]=i;
   };
   #Print the header for the output in Sprime .score-like format:
   print "CHROM", "POS", "ID", "REF", "ALT", "SEGMENT", "ALLELE", "ASSTATE";
}
filenum==2&&FNR>1{
   #Generate the common prefix for all .score file lines at this site:
   prefix=$querycols["CHROM"] OFS $querycols["POS"] OFS $querycols["ID"] OFS $querycols["REF"] OFS $querycols["ALT"];
   for (i in gtcols) {
      #We need the individual allele IDs for each haplotype:
      ploidy=split($gtcols[i], gt, "[/|]");
      for (h=1; h<=ploidy; h++) {
         #Construct the key for the tractcounts/tractindex hashes:
         hapchrom=i"_"h SUBSEP $querycols["CHROM"];
         if (length(debug) > 0 && debug >= 3) {
            print hapchrom, tractindex[hapchrom], tractcounts[hapchrom] > "/dev/stderr";
         };
         #This likely shouldn't ever be the case, but skip if a haplotype
         # doesn't have any tracts on that chromosome:
         if (hapchrom in tractcounts) {
            #Construct the tract key, moving past tracts until the current
            # position is either contained in the tract or is not in a tract:
            tractkey=hapchrom SUBSEP tractindex[hapchrom];
            if (length(debug) > 0 && debug >= 4) {
               print tractkey, $querycols["POS"], tractstarts[tractkey], tractends[tractkey] > "/dev/stderr";
            };
            while ($querycols["POS"] > tractends[tractkey] && tractindex[hapchrom] < tractcounts[hapchrom]) {
               tractindex[hapchrom]++;
               tractkey=hapchrom SUBSEP tractindex[hapchrom];
               if (length(debug) > 0 && debug >= 4) {
                  print tractkey, $querycols["POS"], tractstarts[tractkey], tractends[tractkey] > "/dev/stderr";
               };
            };
            #If the current position is in a tract, output a .score line
            # for the current haplotype:
            if ($querycols["POS"] >= tractstarts[tractkey] && $querycols["POS"] <= tractends[tractkey]) {
               print prefix, tractids[tractkey], gt[h], ASstates[tractkey];
            };
         };
      };
   };
}
