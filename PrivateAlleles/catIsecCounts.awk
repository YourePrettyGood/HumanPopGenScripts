#!/bin/awk -f
#This script takes the output of multiple isecAlleleCounts.awk runs
# and combines them, outputting columns for all possible combinations
# of input files for the intersections, as some isecAlleleCounts.awk
# outputs might be missing columns if those columns are 0.
#Required arguments:
#  numfiles: The number of files input into the isec call
function tobitstring(x,n,origx,bs,i,bit) {
   origx=x;
   bs="";
   for (i=n-1; i>=0; i--) {
      bit=2^i;
      if (x - bit >= 0) {
         bs=bs"1";
         x-=bit;
      } else {
         bs=bs"0";
      };
   };
   if (x != 0) {
      print "Error converting "origx" to bitstring, remainder "x" after "n" iterations getting "bs"." > "/dev/stderr";
   };
   return bs;
}
BEGIN{
   FS="\t";
   OFS=FS;
   #Check for required argument:
   if (length(numfiles) == 0 || numfiles <= 0) {
      print "Missing or non-positive numfiles, please set it appropriately. Quitting." > "/dev/stderr";
      exit 2;
   };
   #Set up a hash of all the possible file combinations:
   maxbs=2^numfiles - 1;
   if (length(debug) > 0 && debug >= 1) {
      print "numfiles="numfiles, "maxbs="maxbs > "/dev/stderr";
   };
   for (i=0; i<=maxbs; i++) {
      groupcols[tobitstring(i,numfiles)]=i;
      if (length(debug) > 0 && debug >= 2) {
         print "groupcols["tobitstring(i, numfiles)"] set to "i > "/dev/stderr";
      };
   };
   #Make sure we iterate over the file combinations in the same order:
   PROCINFO["sorted_in"]="@val_num_desc";
   #Print out a header:
   printf "CHROM\tSubset\tVariantType";
   for (i in groupcols) {
      printf "\t%s", i;
   };
   printf "\n";
}
FNR==1{
   #Don't let previous files contaminate the current file:
   delete colmap;
   #Generate a map from file combinations to extant columns:
   for (i=1; i<=NF; i++) {
      if ($i in groupcols) {
         colmap[$i]=i;
         if (length(debug) > 0 && debug >= 2) {
            print FILENAME" colmap["$i"] set to "i > "/dev/stderr";
         };
      };
   };
}
FNR>1{
   #Print the fixed columns:
   printf "%s\t%s\t%s", $1, $2, $3;
   #Now print the file combination columns:
   for (i in groupcols) {
      #If the combination exists, print it's allele count, else print 0:
      if (i in colmap) {
         printf "\t%i", $colmap[i];
      } else {
         printf "\t0";
      };
   };
   printf "\n";
}
