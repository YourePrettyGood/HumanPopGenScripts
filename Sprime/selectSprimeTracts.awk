#!/bin/awk -f
#This script extracts lines from a file based on matching of
# keys to the values in a key column.
#The input file should be a tab-separated file with a header
# row, and the value of keycol must be found in that header.
#The original purpose of this script was to extract Sprime
# tracts for a particular target population from a summary
# BED file with header.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(keycol) == 0) {
      print "Missing keycol, please specify it. Quitting.\n" > "/dev/stderr";
      exit 2;
   };
   if (length(key) == 0) {
      print "Missing key, assuming you want all input tracts.\n" > "/dev/stderr";
   };
   #Construct a hash to quickly search if any keys in the input list
   # are found:
   n_keys=split(key, keyarr, ",");
   for (i=1; i<=n_keys; i++) {
      keys[keyarr[i]]=i;
   };
}
#Keep track of the column names from the header line:
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Make sure the key column exists:
   if (!(keycol in cols)) {
      print "Did not find keycol "keycol" in header/colnames. Please check the input file. Quitting.\n" > "/dev/stderr";
      exit 3;
   };
   print;
}
#Only output tracts if their key is in the input list of keys:
NR>1{
   #If no key list was specified, just output all input tracts:
   if (length(key) == 0) {
      print;
   } else if ($cols[keycol] in keys) {
      print;
   };
}
