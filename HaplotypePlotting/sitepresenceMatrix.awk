#!/bin/awk -f
#This script takes an rsid to hs37d5 coordinate conversion table for
# the previously-published SNPs, lists of the previously-published SNPs
# (i.e. Banday et al. 2022 Nature Genetics, Gittelman et al. 2016 Cell,
# and Vespasiani et al. 2022 PLoS Genetics), a list of emVars from the
# MPRA, and lists of core haplotype SNPs and Sprime sites, followed
# by a final file of the site coordinates to consider for the output
# matrix.
#The output is a binary matrix indicating whether or not a given site
# was present in a given category/input file.
#With the `no_published` option set to 1, only the emVar list, a single
# core haplotype site list, an S' site list, and the final file of
# site coordinates to consider for the output matrix.
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of which file we're on:
   filenum=0;
   #Print a header line:
   if (length(no_published) > 0 && no_published == "1") {
      print "CHROM:POS", "Core", "Sprime", "emVar";
   } else {
      print "CHROM:POS", "CoreNea", "CoreDen", "SprimeDen", "Banday", "Gittelman", "Vespasiani", "emVar";
   };
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
filenum==1{
   if (length(no_published) > 0 && no_published == "1") {
      emvar[$0]=1;
   } else {
      rs[$1]=$2 OFS $3;
   };
}
filenum==2{
   if (length(no_published) > 0 && no_published == "1") {
      core[$0]=1;
   } else {
      banday[rs[$1]]=1;
   };
}
filenum==3{
   if (length(no_published) > 0 && no_published == "1") {
      sprime[$0]=1;
   } else {
      gittelman[rs[$1]]=1;
   };
}
filenum==4{
   if (length(no_published) > 0 && no_published == "1") {
      c=$0 in core ? "1" : "0";
      s=$0 in sprime ? "1" : "0";
      e=$0 in emvar ? "1" : "0";
      print $1":"$2, c, s, e;
   } else {
      vesp[rs[$1]]=1;
   };
}
filenum==5{
   emvar[$0]=1;
}
filenum==6{
   corenea[$0]=1;
}
filenum==7{
   coreden[$0]=1;
}
filenum==8{
   sprime[$0]=1;
}
filenum==9{
   cn=$0 in corenea ? "1" : "0";
   cd=$0 in coreden ? "1" : "0";
   b=$0 in banday ? "1" : "0";
   g=$0 in gittelman ? "1" : "0";
   v=$0 in vesp ? "1" : "0";
   e=$0 in emvar ? "1" : "0";
   sd=$0 in sprime ? "1" : "0";
   print $1":"$2, cn, cd, sd, b, g, v, e;
}
