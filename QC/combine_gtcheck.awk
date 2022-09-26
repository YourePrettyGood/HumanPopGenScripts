#!/bin/awk -f
#This script combines the outputs from per-chromosome runs of
# bcftools gtcheck, so we get a genome-wide estimate of the
# pairwise genotype discordance between each sample in the
# query and each sample in the ground truth ("target").
#We calculated for every pair, but the inputs need not be
# for all pairs, so the output is dependent on the particular
# bcftools gtcheck command used for the inputs.
#The output is in ascending order of the number of mismatches,
# so you might want to do an extra sort -k5,5g on the
# output to sort in ascending order of genotype discordance.
#Such a sort shouldn't affect much, as the number of
# compared sites shouldn't vary too much between samples.
BEGIN{
   FS="\t";
   OFS=FS;
}
#Sum up the numbers of mismatches (mm) and the numbers of
# compared sites (c) for each query-target sample pair:
/^DC/{
   mm[$2,$3]+=$4;
   c[$2,$3]+=$6;
}
END{
   #Output in descending order of mismatches/discordance,
   # so highest concordance first:
   PROCINFO["sorted_in"]="@val_num_desc";
   #Output for each query-target sample pair:
   #Last column is the genotype discordance rate for that individual.
   for (i in mm) {
      split(i, a, SUBSEP);
      print a[1], a[2], mm[i], c[i], mm[i]/c[i];
   };
}
