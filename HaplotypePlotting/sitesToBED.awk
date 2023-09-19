#!/bin/awk -f
#Simply takes a list of sites (two tab-separated columns: CHROM and POS, no header)
# and converts it into BED format.
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   print $1, $2-1, $2;
}
