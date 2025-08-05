#!/bin/awk -f
#This script takes the STDOUT of yak trioeval and outputs two files listing
# the contigs names in an assembly that should be kept or swapped based on
# trio paternal kmer matching. The concept is pretty simple: Calculate
# the proportion of parental kmer matches corresponding to the opposite of
# the target haplotype and compare this proportion to a fixed threshold
# (default = 0.7). If a contig has too high of a proportion of kmer matches
# to the wrong parent, it should get swapped, otherwise it should be kept.
#The two output files are named with _hap[12]_ctgs.txt suffixes and are
# appended to so that you can run this script sequentially on the hap1 and
# hap2 trioeval runs to produce two lists of contigs for the final assemblies.
#Options:
# outprefix:  Prefix for the two output files (required)
# swapthresh: Threshold of the proportion of wrong parental kmers for swapping
#             (default: 0.7)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(outprefix) == 0) {
      print "No outprefix supplied, please provide a value for outprefix. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(swapthresh) == 0) {
      swapthresh=0.7;
      print "No threshold for swapping provided, using 0.7." > "/dev/stderr";
   };
   if (swapthresh < 0.0 || swapthresh > 1.0) {
      print "Invalid threshold for swapping provided, must be between 0 and 1. Quitting." > "/dev/stderr";
      exit 3;
   };
}
/^S/{
   denom=$3+$4;
   if (denom > 0) {
      matratio[$2]=$3/denom;
   } else {
      #If no kmer info, default to 0 so neither direction indicates a swap:
      matratio[$2]=0.0;
   };
}
/^N/{
   if ($2 > $3) {
      #More maternal kmers than paternal, so filter on patratio=1-matratio:
      inferredhap="1";
      keep_ctgs=outprefix"_hap1_ctgs.txt";
      swap_ctgs=outprefix"_hap2_ctgs.txt";
   } else {
      #More paternal kmers than maternal, so filter on matratio:
      #This also catches the rare case where the paternal kmer counts are equal,
      # but isn't necessarily the right way to handle it.
      inferredhap="2";
      keep_ctgs=outprefix"_hap2_ctgs.txt";
      swap_ctgs=outprefix"_hap1_ctgs.txt";
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (ctg in matratio) {
      if (inferredhap == "1") {
         testratio=1.0-matratio[ctg];
      } else {
         testratio=matratio[ctg];
      };
      if (testratio >= swapthresh) {
         print ctg >> swap_ctgs;
      } else {
         print ctg >> keep_ctgs;
      };
   };
}
