#!/bin/awk -f
#This script takes a set of GFF3 lines (or an entire GFF3) and outputs one
# particular tag's value for a particular feature type.
#This is generally useful for things like extracting gene IDs from the
# output of overlappingFeatures.pl, or equivalently extracting transcript
# IDs or even exon IDs if available, gene names, etc.
#Optional arguments:
# tag:     GFF3 tag name to extract
#          (default: ID)
# feature: GFF3 feature type to examine (column 3 of GFF3)
#          (default: gene)
# trim:    Flag to trigger trimming of "[.]\d+$" from the end of the tag value
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(tag) == 0) {
      tag="ID";
   };
   if (length(feature) == 0) {
      feature="gene";
   };
   skipFASTA=0;
}
/^##FASTA/{
   skipFASTA=1;
}
!/^#/&&!skipFASTA{
   if ($3 != feature) {
      next;
   };
   split($9, tags, "[ ]?;[ ]?");
   for (i in tags) {
      split(tags[i], tagelems, "=");
      if (tagelems[1] == tag) {
         if (length(trim) > 0) {
            sub("[.][0-9]+$", "", tagelems[2]);
         };
         print tagelems[2];
      };
   };
}
