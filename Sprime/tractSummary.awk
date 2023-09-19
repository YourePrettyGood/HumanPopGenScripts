#!/bin/awk -f
#This script takes a BED file and calculates the sum of
# interval lengths as well as the count of intervals
# and outputs the sum, count, and average interval length.
#This output is prefixed with two input variables:
# origin (the tract origin)
# pop or region (the target population or region)
#These are custom to the initial use case, which is
# summarizing reconstructed archaic sequence from
# Sprime tract projections.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(origin) == 0) {
      print "Missing origin variable, please specify it. Quitting.\n" > "/dev/stderr";
      exit 2;
   };
   if (length(pop) == 0 && length(region) == 0) {
      print "Missing pop or region variable, please specify one. Quitting.\n" > "/dev/stderr";
      exit 3;
   };
   if (length(pop) > 0 && length(region) > 0) {
      print "Both pop and region are set, please only use one. Quitting.\n" > "/dev/stderr";
      exit 4;
   };
   #By default, we assume that we're being fed parts of a tiling path:
   #This only affects output header contents.
   #Another possible type would be "Projection" when projections
   # haven't been merged.
   if (length(type) == 0) {
      type="Path";
   };
   #Set initial values so we don't ever print out blank strings:
   sum=0;
   count=0;
   #If a header is requested, output it:
   if (length(header) > 0 && header > 0) {
      if (length(pop) > 0) {
         print "TractOrigin", "Population", "Total"type"Length", "Num"type"Parts", "AveragePartLength";
      } else if (length(region) > 0) {
         print "TractOrigin", "Region", "Total"type"Length", "Num"type"Parts", "AveragePartLength";
      };
   };
}
!/^#/{
   sum+=$3-$2;
   count+=1;
}
END{
   if (count == 0) {
      avglen="NA";
   } else {
      avglen=sum/count;
   };
   if (length(pop) > 0) {
      print origin, pop, sum, count, avglen;
   } else if (length(region) > 0) {
      print origin, region, sum, count, avglen;
   };
}
