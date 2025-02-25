#!/bin/awk -f
#This script extracts windows deemed "significant" from the U and Q95
# statistics output from the UQ95.nf pipeline for selected group B
# using outgroup A. Windows are deemed "significant" if their U and
# Q95 values are BOTH greater than values determined from the first
# input file.
#These thresholds are determined on a per-archaic-origin basis, and
# the first input file may be determined by e.g. the 95th percentile
# of genome-wide U and Q95 values for fixed threshold frequencies
# w and x.
#The first input file can be created with the following R code:
# library(tidyverse)
# #Set the quantile for genome-wide significance:
# gwthresh <- 0.95
# #Evaluate the genome-wide thresholds for U and Q95:
# UQ95 <- read_tsv('[path to output of UQ95.nf pipeline]')
# UQ95 %>% group_by(Origin, A, B) %>%
#    summarize(Uthresh=quantile(U_ABCD, probs=gwthresh, na.rm=TRUE),
#              Q95thresh=quantile(Q95_ABCD, probs=gwthresh, na.rm=TRUE)) %>%
#    write_tsv('[path to thresholds output file]', quote="none", escape="none")
#
#For details on the parameters of U and Q95 (i.e. A, B, C, D, w, x, y, z),
# see Racimo et al. 2017 MBE.
#Arguments:
# A:         Outgroup population(s), comma-separated string if multiple.
#            Period (.) can be used as a wildcard. A is exactly matched
#            to the input file contents.
# B:         Target population(s), comma-separated string if multiple.
#            Period (.) can be used as a wildcard. B is exactly matched
#            to the input file contents.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(A) == 0) {
      print "Missing A argument, please set it." > "/dev/stderr";
      print "A is the outgroup of the U and Q95 statistics to extract." > "/dev/stderr";
      exit 2;
   };
   if (length(B) == 0) {
      print "Missing B argument, please set it." > "/dev/stderr";
      print "B is the target group of the U and Q95 statistics to extract." > "/dev/stderr";
      exit 3;
   };
   filter=1;
   filenum=0;
}
FNR==1{
   filenum+=1;
}
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      threshcols[$i]=i;
   };
}
filenum==1&&FNR>1{
   if ($threshcols["A"] == A && $threshcols["B"] == B) {
      Uthresh[$threshcols["Origin"]]=$threshcols["Uthresh"];
      Q95thresh[$threshcols["Origin"]]=$threshcols["Q95thresh"];
   };
}
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum==2&&FNR>1{
   if ($cols["A"] == A && $cols["B"] == B) {
      if (filter > 0) {
         if ($cols["U_ABCD"] != "NA" && $cols["Q95_ABCD"] != "NA") {
            if ($cols["U_ABCD"] > Uthresh[$cols["Origin"]] && $cols["Q95_ABCD"] > Q95thresh[$cols["Origin"]]) {
               print $cols["CHROM"], $cols["start"]-1, $cols["end"], $cols["A"], $cols["B"], $cols["Origin"], $cols["U_ABCD"], $cols["Q95_ABCD"];
            };
         };
      } else {
         print $cols["CHROM"], $cols["start"]-1, $cols["end"], $cols["A"], $cols["B"], $cols["Origin"], $cols["U_ABCD"], $cols["Q95_ABCD"];
      };
   };
}
