#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   n_alts=split($4, alleles, ",");
   alleles[0]=$3;
   for(i=0; i<=n_alts; i++) {
      AC[i]=0;
   };
   for (i=5; i<=NF; i++) {
      ploidy=split($i, gt, "[/|]");
      for (j=1; j<=ploidy; j++) {
         AC[gt[j]]+=1;
      };
   };
   if (AC[0]*AC[1] > 0) {
      print $1, $2-1, $2;
   };
}
