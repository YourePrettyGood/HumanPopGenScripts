#!/bin/awk -f
#This script is a bit of a hacky way to make an hg19 BAM header compatible
# with hs37d5 (at least as far as the major chromosomes are concerned).
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1{
   scafmap[$2]=$1;
}
filenum==2&&/^@/&&!/^@SQ/{
   print;
}
filenum==2&&/^@SQ/{
   for (i=2; i<=NF; i++) {
      n=split($i, snbits, ":");
      if (snbits[1] == "SN" && snbits[2] in scafmap) {
         $i="SN:"scafmap[snbits[2]];
         break;
      };
   };
   print $0;
}
filenum==2&&!/^@/{
   if ($3 in scafmap) {
      $3=scafmap[$3];
   };
   if ($7 in scafmap) {
      $7=scafmap[$7];
   };
   print $0;
}
