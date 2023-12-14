#!/bin/awk -f
#A simple script to extract sample IDs in the offspring->father->mother
# order expected by generate_multihetsep.py with --trio 0,1,2.
#The input is a ped file with format as indicated by
# the first six columns of:
# https://github.com/samtools/bcftools/blob/44deedcd1a56506191e95ebb35a549b0b6dbc9a1/plugins/trio-stats.c#L183
#We basically just filter on the first column (familyID) and output
# sampleID, then paternalID, then maternalID.
#Required arguments:
# trios: A comma-separated list of trio familyIDs to extract
#Optional arguments:
# no_O:  Whether to exclude the offspring from the output
#        (default: Include offspring)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(trios) == 0) {
      print "trios list is empty, cannot proceed. Quitting." > "/dev/stderr";
      exit 2;
   };
   n_trios=split(trios, trioarr, ",");
   for (i=1; i<=n_trios; i++) {
      keeptrios[trioarr[i]]=i;
   };
}
#Make sure we skip the header if it exists:
!/^familyID/{
   if ($1 in keeptrios) {
      if (length(no_O) == 0 || no_O < 1) {
         print $2;
      };
      print $3;
      print $4;
   };
}
