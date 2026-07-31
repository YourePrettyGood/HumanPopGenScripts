#!/bin/awk -f
#This script takes a BEDGRAPH file of mask depth and outputs
# intervals with depth less than or equal to a threshold
# defined by maxfrac * n.
#This comparison can be inverted by specifying negate.
#Mandatory arguments:
# maxfrac: Maximum fraction of samples masked to call as passing
#          (Must be in [0,1])
# n:       Total number of samples
#Optional arguments:
# negate:  Only output intervals with depth > threshold
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(maxfrac) == 0) {
      print "Missing maxfrac, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(n) == 0) {
      print "Missing n, please specify it." > "/dev/stderr";
      exit 3;
   };
   if (maxfrac < 0 || maxfrac > 1) {
      print "Invalid maxfrac, must satisfy 0 <= maxfrac <= 1." > "/dev/stderr";
      exit 4;
   };
   if (n < 1) {
      print "Invalid n, must satisfy n > 0." > "/dev/stderr";
      exit 5;
   };
   #Make sure the thresholds get treated as numerics:
   maxn=maxfrac*(n+0);
   #Make sure negate takes on a usable value:
   #Anything other than absent, 0, or variations on No or False is treated
   # as true.
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "", negate);
   sub(/^[Nn][Oo]?$/, "", negate);
   sub(/^0$/, "", negate);
   sub(/^[Tt]([Rr][Uu][Ee])?$/, "1", negate);
   sub(/^[Yy]([Ee][Ss])?$/, "1", negate);
   if (length(negate) == 0) {
      negate=0;
   };
}
{
   if ((!negate && $4 <= maxn) || (negate && $4 > maxn)) {
      print $1, $2, $3;
   };
}
