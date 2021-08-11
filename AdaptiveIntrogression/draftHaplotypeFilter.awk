#!/bin/awk -f
#This script takes the output of r2toDraftHaplotypes.awk and filters
# on both the number of SNPs in the cluster and the fraction of
# missing edges.
#The idea is to replicate the first clustering step of Gittelman et al.
# 2016 Current Biology to generate "initial" (draft) haplotypes for
# estimating introgressed haplotype frequencies. In particular, this
# deals with both the explicit filter of minimum 3 tag SNPs, and an
# implicit filter that the haplotype should be sufficiently linked
# across the entire haplotype.
#Optional arguments:
# nsnps:        The minimum number of SNPs in a cluster (default: 3)
# missingness:  The maximum fraction of edges missing for a cluster
#               (default: 0.0)
#Note that setting missingness too low will tend to filter out the
# longest haplotypes, as LD decay will almost guarantee some edges
# fall under the r^2 threshold.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(nsnps) == 0) {
      nsnps=3;
   };
   if (length(missingness) == 0) {
      missingness=0.0;
   };
}
#Retain the header, and add the missingness column:
NR==1{
   print $0, "MISSINGNESS";
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
#Filter the clusters:
NR>1{
   #The total number of edges in a complete graph of n nodes is n*(n-1)/2,
   # so the missing fraction is just the number of missing edges over this:
   missingfraction=$cols["MISSINGEDGES"]/($cols["CLUSTERSIZE"]*($cols["CLUSTERSIZE"]-1)/2);
   if ($cols["CLUSTERSIZE"] >= nsnps && missingfraction <= missingness) {
      print $0, missingfraction;
   };
}
