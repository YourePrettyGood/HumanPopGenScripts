#!/bin/awk -f
#This script does simple filtering of the table of Sprime tract frequencies
# by both population and minimum tract frequency.
#We assume that the input has a header line with column names, and there are
# columns as follows:
# Population = name/ID of the population run through Sprime
# Chromosome = chromosome where the tract was identified
# BEDStart =   0-based position of the start of the tract (i.e. the left-most
#              position in the tract)
# End =        1-based position of the end of the tract (i.e. the right-most
#              position in the tract)
# AAF =        Tract frequency (AAF stands for "Archaic Allele Frequency",
#              as I've been calculating tract frequency by taking the median
#              of the frequencies of the alleles in the tract)
#Required options:
# pop:    ID of the population whose tracts you want to filter
# thresh: Minimum tract frequency to retain/emit (range: float from 0 to 1)
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   if ($cols["Population"] == pop && $cols["AAF"] >= thresh) {
      print "chr"$cols["Chromosome"], $cols["BEDStart"], $cols["End"];
   };
}
