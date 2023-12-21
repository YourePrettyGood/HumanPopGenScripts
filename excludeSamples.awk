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
      samplecol=1;
      if (length(samplecolname) > 0) {
         for (i=1; i<=NF; i++) {
            if ($i == samplecolname) {
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
