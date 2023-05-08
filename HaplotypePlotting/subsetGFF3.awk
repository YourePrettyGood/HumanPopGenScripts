#!/bin/awk -f
#This script takes a GFF3 as input and outputs header lines and only a
# subset of records based on the value of the feature column (column 3).
#It also optionally trims the "chr" prefix from the scaffold column
# (column 1).
#It's kind of a silly workaround for using the GENCODE annotation with
# hs37d5.
#
#Optional arguments:
# feature: GFF3 feature type to retain (column 3 of GFF3)
#          (default: exon)
# trimchr: Flag to trigger trimming the "chr" prefix from column 1
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(feature) == 0) {
      feature="exon";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", trimchr);
   sub(/^[Nn][Oo]?$/, "0", trimchr);
}
#Feed all header lines through:
/^#/{
   print;
}
#Only retain record lines if their feature matches:
!/^#/&&$3==feature{
   if (length(trimchr) > 0 && trimchr > 0) {
      sub("chr", "", $1);
   };
   print $0;
}
