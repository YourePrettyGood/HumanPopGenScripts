#!/bin/awk -f
#This script takes a bunch of individual .sample files and merges them
# for use with Relate. The general idea is to merge if you have
# per-chromosome .sample files from converting VCF to .haps/.sample.
#Relate ignores the missing data proportion column, so we don't
# do anything special to combine these values across files.
#We do check that the individuals match across files.
BEGIN{
   OFS=" ";
   n_samples=0;
   #Error state to fail out at the end if there was a mismatch:
   errstate=0;
   #Keep track of which file we're on:
   filenum=0;
}
FNR==1{
   #Keep track of which file we're on:
   filenum++;
}
#Pass through the header:
filenum==1&&FNR==1{
   print;
}
#Pass through the first file's sample lines and store:
filenum==1&&FNR>1{
   n_samples+=1;
   found[$1,$2]=n_samples;
   print;
}
#Check other files for matching sample lines:
filenum>1&&FNR>1{
   if (!(($1,$2) in found)) {
      #n_samples+=1;
      print "Mismatched sample line in file "filenum" ("FILENAME"): "$0 > "/dev/stderr";
      errstate+=1;
   };
}
END{
   if (errstate > 0) {
      exit 2;
   };
}
