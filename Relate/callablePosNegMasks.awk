#!/bin/awk -f
#This script takes a per-base depth BED and outputs two gzipped
# BED files of the positive and negative depth masks.
#The positive depth mask [sample ID].callable.bed.gz indicates
# the genomic regions where mindp <= DP < supdp.
#The negative depth mask [sample ID].mask.bed.gz indicates the
# genomic regions where DP < mindp or DP >= supdp.
#Mandatory arguments:
# id:    Sample ID
# mindp: Minimum depth to be considered callable
# supdp: Supremum depth to be considered callable
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(id) == 0) {
      print "Missing id, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(mindp) == 0) {
      print "Missing mindp, please specify it." > "/dev/stderr";
      exit 3;
   };
   if (length(supdp) == 0) {
      print "Missing supdp, please specify it." > "/dev/stderr";
      exit 4;
   };
   if (mindp < 0 || supdp <= mindp) {
      print "Invalid mindp and/or supdp, must satisfy 0 <= mindp < supdp." > "/dev/stderr";
      exit 5;
   };
   #Make sure the thresholds get treated as numerics:
   mindp=mindp+0;
   supdp=supdp+0;
   #Set up the output pipes:
   callable="bedtools merge -i - | gzip -9 > "id".callable.bed.gz";
   mask="bedtools merge -i - | gzip -9 > "id".mask.bed.gz";
}
{
   if ($4 >= mindp && $4 < supdp) {
      #Callable site:
      print $1, $2, $3 | callable;
   } else {
      print $1, $2, $3 | mask;
   };
}
