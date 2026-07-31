#!/bin/awk -f
#This script simply takes the asymmetric complement of the samples in
# the second input file relative to those of the first. That is, any
# samples in the first file are excluded from the second file (if found).
#I've also added the capacity to handle a more complex second file
# consisting of a header line and multiple columns.
#The header is retained and used to identify the sample ID column,
# and sample lines are retained if their ID is not in the exclude list.
BEGIN{
   FS="\t";
   #Note: This doesn't actually convert input spaces to output tabs if
   # spaces are auto-detected, since we use print; instead of print $0;
   # below.
   OFS=FS;
   if (length(header) == 0) {
      header=0;
   };
   if (length(samplecolname) == 0 && length(samplecol) == 0) {
      samplecol=1;
   };
   if (length(negate) == 0) {
      negate=0;
   };
   filenum=0;
}
FNR==1{
   filenum++;
   #Auto-detect space-separation vs. tab-separation:
   #This approach has a bug if the second file has a single column
   # that contains spaces.
   if (filenum==2) {
      test_n_cols=split($1, test_cols, " ");
      if (NF == 1 && test_n_cols > 1) {
         FS=" ";
      };
   };
}
filenum==1{
   exclude[$1]=1;
}
filenum==2{
   if (FNR==1 && header > 0) {
#      samplecol=1;
      if (length(samplecolname) > 0) {
         #Per https://ftp.gnu.org/old-gnu/Manuals/gawk-3.1.1/html_node/Field-Splitting-Summary.html#Field%20Splitting%20Summary
         # changing FS does not cause a reparse of the same line for
         # POSIX-compliant awk implementations, so we need to do the
         # reparse explicitly ourselves:
         reparse_NF=split($0, reparse_cols, FS);
         for (i=1; i<=reparse_NF; i++) {
            if (reparse_cols[i] == samplecolname) {
               samplecol=i;
            };
         };
      };
      print;
   } else if (negate > 0) {
      if ($samplecol in exclude) {
         print;
      };
   } else if (!($samplecol in exclude)) {
      print;
   };
}
