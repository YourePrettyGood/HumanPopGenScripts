#!/bin/awk -f
#This script extracts windows deemed "significant" from the U and Q95
# statistics output from the UQ95.nf pipeline for selected group B
# using outgroup A. Windows are deemed "significant" if their U value
# is greater than the Uthresh argument AND if their Q95 value is greater
# than the Q95thresh argument.
#Arguments:
# A:         Outgroup population(s), comma-separated string if multiple.
#            Period (.) can be used as a wildcard. A is exactly matched
#            to the input file contents.
# B:         Target population(s), comma-separated string if multiple.
#            Period (.) can be used as a wildcard. B is exactly matched
#            to the input file contents.
# Uthresh:   Threshold U_{A,B,C,D} value above which a window is significant
# Q95thresh: Threshold Q95_{A,B,C,D} value above which a window is significant
#For details on the parameters of the U and Q95 statistics (i.e. A, B, C, D, w,
# x, y, z), see Racimo et al. 2017 MBE
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(A) == 0) {
      print "Missing A argument, please set it." > "/dev/stderr";
      print "A is the outgroup of the U and Q95 statistics to extract." > "/dev/stderr";
      exit 2;
   };
   if (length(B) == 0) {
      print "Missing B argument, please set it." > "/dev/stderr";
      print "B is the target group of the U and Q95 statistics to extract." > "/dev/stderr";
      exit 3;
   };
   filter=1;
   if (length(Uthresh) == 0 && length(Q95thresh) == 0) {
      filter=0;
   };
   if (length(Uthresh) == 0) {
      Uthresh=-1;
   };
   if (length(Q95thresh) == 0) {
      Q95thresh=-1;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   if ($cols["A"] == A && $cols["B"] == B) {
      if (filter > 0) {
         if ($cols["U_ABCD"] != "NA" && $cols["Q95_ABCD"] != "NA") {
            if ($cols["U_ABCD"] > Uthresh && $cols["Q95_ABCD"] > Q95thresh) {
               print $cols["CHROM"], $cols["start"]-1, $cols["end"], $cols["A"], $cols["B"], $cols["Origin"], $cols["U_ABCD"], $cols["Q95_ABCD"];
            };
         };
      } else {
         print $cols["CHROM"], $cols["start"]-1, $cols["end"], $cols["A"], $cols["B"], $cols["Origin"], $cols["U_ABCD"], $cols["Q95_ABCD"];
      };
   };
}
