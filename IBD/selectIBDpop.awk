#!/bin/awk -f
#This script retains only IBD (or HBD) tracts in the Browning .ibd format
# (e.g. produced by Refined-IBD or hap-IBD) between samples identified
# in the first input file.
#Thus, the inputs are:
# File 1) A list of sample IDs, one per line, to retain
# File 2) The uncompressed .ibd (or .hbd) file
#The output is a filtered version of the second input.
#Options:
# ind1col: Column number for the ID of individual 1 of the IBD tract
#          (Default: 1)
# ind2col: Column number for the ID of individual 2 of the IBD tract
#          (Default: 3)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(ind1col) == 0) {
      ind1col=1;
      print "Defaulting ind1col to column 1" > "/dev/stderr";
   };
   if (length(ind2col) == 0) {
      ind2col=3;
      print "Defaulting ind2col to column 3" > "/dev/stderr";
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
   keep[$1]=1;
}
filenum==2{
   if ($ind1col in keep && $ind2col in keep) {
         print;
   };
}
