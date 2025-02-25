#!/bin/awk -f
#This script takes in the *decoded.*.txt outputs from hmmix decode and
# converts any tracts in the Archaic HMM state with sufficiently high
# posterior probability (i.e. >= minposterior) into BED format.
#This is useful for computing overlaps with introgressed tract inferences
# from other methods (e.g. ArchaicSeeker2.0, Sprime projections, S*)
# using something like BEDtools.
#Arguments:
# minposterior: Minimum threshold for posterior probability of an Archaic tract
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(minposterior) == 0) {
      print "minposterior variable is missing, please set it." > "/dev/stderr";
      print "minposterior is the minimum threshold of the posterior probability of HMM state for an Archaic tract." > "/dev/stderr";
      exit 2;
   };
   if (minposterior < 0.0 || minposterior > 1.0) {
      print "minposterior "minposterior" is outside of [0,1], please fix it." > "/dev/stderr";
      exit 3;
   };
}
#Map column names to indices:
FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
FNR>1{
   if ($cols["state"] == "Archaic" && $cols["mean_prob"] >= minposterior) {
      print $cols["chrom"], $cols["start"], $cols["end"];
   };
}
