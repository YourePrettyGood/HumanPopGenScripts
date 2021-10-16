#!/bin/awk -f
#This script takes the output of bedtools intersect -wao with -a as a BED6
# of S' tracts (TractID in column 4 at a minimum) and -b as a GFF3 (or GFF3
# lines) and outputs a TSV of the tract ID, chromosome, start, end, and a
# comma-separated list of (by default) gene IDs overlapping that tract.
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
   #Expected number of columns from the -a file:
   #CHROM, START, END, NAME, SCORE, STRAND
   FNF=6
   #Expected number of columns from the -b file:
   #CHROM, SOURCE, FEATURE, START, END, SCORE, STRAND, PHASE, ATTR
   SNF=9
}
{
   #Only examine overlaps with the desired feature type:
   if ($(FNF+3) != feature) {
      next;
   };
   #Extract the tract ID from the BED line:
   tractid="";
   split($4, tags, "[ ]?;[ ]?");
   for (i in tags) {
      split(tags[i], tagelems, "=");
      if (tagelems[1] == "TractID") {
         tractid=tagelems[2];
         break;
      };
   };
   tractcoords[tractid]=$1 OFS $2 OFS $3;
   #Extract the feature ID from the GFF3 line, and trim the suffix if requested:
   split($(FNF+9), tags, "[ ]?;[ ]?");
   for (i in tags) {
      split(tags[i], tagelems, "=");
      if (tagelems[1] == tag) {
         if (length(trim) > 0) {
            sub("[.][0-9]+(_[0-9]+)?$", "", tagelems[2]);
         };
         if (tractid in genelist) {
            genelist[tractid]=genelist[tractid]","tagelems[2];
         } else {
            genelist[tractid]=tagelems[2];
         };
      };
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (tractid in genelist) {
      print tractcoords[tractid], tractid, genelist[tractid];
   };
}
