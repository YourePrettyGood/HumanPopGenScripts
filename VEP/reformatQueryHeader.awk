#!/bin/awk -f
#This script takes the output of bcftools query -H or bcftools +split-vep -H and
# reformats the header into a format better for import into R.
BEGIN{
   FS="\t";
   OFS=FS;
}
#Have to put the block for the rest of the lines first to avoid double-printing the header:
!/^#/{
   print;
}
#Reformat the header:
/^#/{
   #Remove the [#] prefixes of each column name:
   gsub("[[0-9]+]", "");
   #Get rid of the # prefix for the header line:
   $0=substr($0, 2);
   #Print it:
   print $0;
}
