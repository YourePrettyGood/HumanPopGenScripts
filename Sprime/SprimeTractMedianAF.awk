#!/bin/awk -f
#This script takes the output of SprimeArchaicAFs.awk and summarizes the
# allele frequency of each tract for each query population using the median,
# minimum, and maximum per-site Sprime-identified archaic allele frequency.
#The minimum and maximum are mainly for diagnostic purposes.
#Optional arguments:
# debug:   Output debug info to STDERR
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   #For each tract, keep track of the interval (start=min(POS), end=max(POS)):
   if ($3 in tractchrom) {
      if ($2 < tractstart[$3]) {
         tractstart[$3]=$2;
      };
      if ($2 > tractend[$3]) {
         tractend[$3]=$2;
      };
   } else {
      tractchrom[$3]=$1;
      tractstart[$3]=$2;
      tractend[$3]=$2;
   };
   #Only store allele frequencies if they aren't NA:
   if ($6 != "NA") {
      #Store them in a comma-separated list:
      if (($7,$3) in tractafs) {
         tractafs[$7,$3]=tractafs[$7,$3]","$6;
      } else {
         tractafs[$7,$3]=$6;
      };
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (i in tractafs) {
      split(i, qpoptract, SUBSEP);
      n=split(tractafs[i], aflist, ",");
      asort(aflist, afs, "@val_num_asc");
      if (length(debug) > 0) {
         printf "%s %s\t%i items", qpoptract[1], qpoptract[2], n > "/dev/stderr";
         for (j=1; j<=n; j++) {
            printf ",%f", afs[j] > "/dev/stderr";
         };
         printf "\n" > "/dev/stderr";
      };
      if (n%2==1) {
         median=afs[(n+1)/2];
      } else {
         median=(afs[n/2]+afs[(n/2)+1])/2;
      };
      print tractchrom[qpoptract[2]], tractstart[qpoptract[2]]-1, tractend[qpoptract[2]], qpoptract[2], qpoptract[1], median, afs[1], afs[n];
   };
}
