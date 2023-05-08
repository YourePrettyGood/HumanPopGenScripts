#!/bin/awk -f
#This script takes the output of bedtools intersect -wao with -a as a BED3
# of sites (e.g. the sites of a haplotype plot) and -b as a GFF3 (or GFF3
# lines) and outputs a TSV of the sites and the feature (if any) that
# overlaps them, as well as the strand of the feature.
#Keep in mind that you will need to deduplicate the output in the case
# of multiple features overlapping the same site.
#
#Optional arguments:
# tag:     GFF3 tag name to extract
#          (default: gene_name)
# feature: GFF3 feature type to examine (column 3 of GFF3)
#          (default: exon)
# trim:    Flag to trigger trimming of "[.]\d+$" from the end of the tag value
# type:    Which gene_type to filter for
#          (no filter applied if variable not set)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(tag) == 0) {
      tag="gene_name";
   };
   if (length(feature) == 0) {
      feature="exon";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", trim);   
   sub(/^[Nn][Oo]?$/, "0", trim);
   #Expected number of columns from the -a file:
   #CHROM, START, END
   FNF=3
   #Expected number of columns from the -b file:
   #CHROM, SOURCE, FEATURE, START, END, SCORE, STRAND, PHASE, ATTR
   SNF=9
   #Print a header:
   print "CHROM:POS", "Feature", "Strand";
}
{
   #Only examine overlaps with the desired feature type:
   if ($(FNF+3) != feature) {
      next;
   };
   #Extract the feature ID from the GFF3 line, and trim the suffix if requested:
   genename="";
   genetype="";
   split($(FNF+9), tags, "[ ]?;[ ]?");
   for (i in tags) {
      split(tags[i], tagelems, "=");
      if (tagelems[1] == tag) {
         if (length(trim) > 0 && trim > 0) {
            sub("[.][0-9]+(_[0-9]+)?$", "", tagelems[2]);
         };
         genename=tagelems[2];
      };
      if (tagelems[1] == "gene_type") {
         genetype=tagelems[2];
      };
   };
   #Note that checking for genename existence is an indirect
   # way to check if an overlap exists.
   if (genename != "") {
      if (length(type) > 0 && genetype != type) {
         next;
      };
      print $1":"$3, genename, $(FNF+7);
   };
}
