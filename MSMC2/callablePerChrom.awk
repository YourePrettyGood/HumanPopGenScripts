#!/bin/awk -f
#A very simple script to extract the callable intervals for a sample on
# a given chromosome from the .quantized.bed output of mosdepth --quantize.
#This is needed for the --mask argument of generate_multihetsep.py of MSMC
# and MSMC2. generate_multihetsep.py completely ignores the chromosome
# column of the BED, so if you don't split this way, you can get some
# masking artifacts (e.g. all chromosomes getting the chr1 mask intervals
# masked, and even some weirder cases).
#Required arguments:
# chrom: Name of the chromosome to retain
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(chrom) == 0) {
      print "Missing chrom variable, please provide it. Quitting." > "/dev/stderr";
      exit 2;
   };
}
$1==chrom&&$4=="CALLABLE"{
   print $1, $2, $3;
}
