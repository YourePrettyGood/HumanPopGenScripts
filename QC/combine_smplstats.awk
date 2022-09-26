#!/bin/awk -f
#This script combines the results from the per-chromosome runs of
# bcftools +smpl-stats to determine genome-wide values of
# heterozygosity, missingness, and Ts/Tv for each individual.
#The heterozygosity estimate is uncorrected, as it simply
# calculates the number of heterozygotes at all passing sites
# regardless of the number of alleles at that site (so ALT-ALT
# heterozygotes are counted) and divides by the total number
# of passing sites.
BEGIN{
   FS="\t";
   OFS=FS;
}
/^FLT/{
   total[$1,$2]+=$3;
   nonref[$1,$2]+=$4;
   homref[$1,$2]+=$5;
   homalt[$1,$2]+=$6;
   het[$1,$2]+=$7;
   hemi[$1,$2]+=$8;
   missing[$1,$2]+=$12;
   ts[$1,$2]+=$13;
   tv[$1,$2]+=$14;
}
END{
   #Sort in ascending order by the heterozygosity estimate:
   PROCINFO["sorted_in"]="@val_num_asc";
   #Output a header line:
   print "Filter", "Sample", "HomRef", "Het", "HomAlt", "Hemi", "Missing", "NonRef", "Total", "HetRate", "HomRefRate", "HomAltRate", "MissingRate", "Ts", "Tv", "Ts/Tv";
   #Output a line for each sample and filter combination:
   for (i in het) {
      split(i, a, SUBSEP);
      print a[1], a[2], homref[i], het[i], homalt[i], hemi[i], missing[i], nonref[i], total[i], het[i]/total[i], homref[i]/total[i], homalt[i]/total[i], missing[i]/total[i], ts[i], tv[i], ts[i]/tv[i];
   };
}
