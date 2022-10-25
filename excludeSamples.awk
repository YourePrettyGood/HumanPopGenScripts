#!/bin/awk -f
#This script simply takes the asymmetric complement of the samples in
# the second input file relative to those of the first. That is, any
# samples in the first file are excluded from the second file (if found).
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1{
   exclude[$1]=1;
}
filenum==2{
   if (!($1 in exclude)) {
      print;
   };
}
