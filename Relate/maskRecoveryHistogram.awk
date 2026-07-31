#!/bin/awk -f
#This script generates an empirical CDF of the recovery of filtered
# variants at different thresholds for mask depth.
#Input is the output of bedtools intersect -wo -a [filtered variants BED]
# -b [mask depth BEDGRAPH] -g [.genome] -sorted
BEGIN{
   FS="\t";
   OFS=FS;
}
#Only store autosome hits:
$1~/^(chr)?[1-9]/{
   #Keep track of the per-threshold site counts:
   c[$7]+=$8;
   #and the total overlap across all thresholds:
   ov+=$8;
   #as well as the total number of filtered variants:
   t+=$3-$2;
}
END{
   PROCINFO["sorted_in"]="@ind_num_desc";
   #Calculate the cumulative site counts with decreasing threshold:
   cumul=0;
   for (i in c) {
      cumul+=c[i];
      print i, c[i], c[i]/ov, cumul/ov, c[i]/t, cumul/t;
   };
}
