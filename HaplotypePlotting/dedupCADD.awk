#!/bin/awk -f
#This script takes the output of addCADD.awk and deduplicates it,
# taking the maximum score at each site.
#Be sure to run sort -k1,1V on the output, otherwise it will be
# in lexicographical order, rather than chromosomal position.
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   print "CHROM:POS", "CADD";;
}
NR>1{
   if ($cols["CADD"] >= cadd[$cols["CHROM:POS"]]) {
      cadd[$cols["CHROM:POS"]]=$cols["CADD"];
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (site in cadd) {
      print site, cadd[site];
   };
}
