#!/bin/awk -f
#This script rearranges the Q file based on the population map used
# for ADMIXTURE/CLUMPAK. Three input files are required:
# 1) The .fam file from PLINK used as input to ADMIXTURE
# 2) The population map file input to the pipeline
# 3) The Q file to rearrange
BEGIN{
   filenum=0;
}
FNR==1{
   filenum++;
}
#First input file is the .fam file, so construct a map from the
# sample ID to the index of the sample in the .fam file:
filenum==1{
   orig[$2]=FNR;
}
#Second input file is the population map, so construct a map from
# the sample ID (column 1) to the index of the sample in the
# population map:
filenum==2{
   reordered[$1]=FNR;
}
#Third input file is the Q file with rows arranged in the .fam
# order, so store the Q file lines:
filenum==3{
   Q[FNR]=$0;
}
END{
   #Iterate over the population map sample IDs in the order they
   # appear, and output the Q file line corresponding to that
   # sample ID:
   #(i here is the sample ID, so orig[i] is the original index)
   PROCINFO["sorted_in"]="@val_num_asc";
   for (i in reordered) {
      print Q[orig[i]];
   };
}
