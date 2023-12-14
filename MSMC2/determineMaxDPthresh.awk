#!/bin/awk -f
#A very simple script to determine the max DP threshold for callable sites
# based on mosdepth.
#The input is the depth dist output of mosdepth
#Required argument:
# quantile: Quantile of DP distribution to determine threshold
#           (e.g. 0.95 would set dpthresh to the 95th percentile of DP)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(quantile) == 0 || quantile < 0 || quantile > 1) {
      print "Missing or invalid quantile specified. Must be in the range [0,1]. Quitting." > "/dev/stderr";
      exit 2;
   };
   maxdpq=1-quantile;
}
$1=="total"&&$3<=maxdpq{
   dpthresh=$2;
}
END{
   print dpthresh;
}
