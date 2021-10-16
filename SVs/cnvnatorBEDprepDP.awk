#!/bin/awk -f
#This script takes the output of bedtools intersect -wao
# -a queryCNVs.bed -b cnvnatorcalls.bed where
# cnvnatorcalls.bed is generated by cnvnatorAnnotatedBED.awk
# and outputs a table ready for R parsing to do genotype
# classification based on DPNORM and overlap
#Output columns:
# CNVID:   The name of the query CNV (provided by the NAME
#          tag in column 4 of the query BED), should be
#          unique in the query BED
# Sample:  Sample ID from the SNAME tag of column 4 of the
#          cnvnator annotated BED
# CNVtype: Like SVTYPE in VCFs
# DPNORM:  Normalized read depth across the window, 1 means
#          the adjusted average (diploid?) depth
# qCNVlen: Length of the query CNV in bp
# Overlap: Length of the query-target overlap in bp
BEGIN{
   FS="\t";
   OFS=FS;
   print "CNVID", "Sample", "CNVtype", "DPNORM", "qCNVlen", "Overlap";
}
NF==13{
   split($4, info, ";");
   for (i in info) {
      split(info[i], tagelems, "=");
      if (tagelems[1] == "NAME") {
         cnvid=tagelems[2];
      };
   };
   qcnvlen=$3-$2;
   split($10, info, ";");
   for (i in info) {
      split(info[i], tagelems, "=");
      if (tagelems[1] == "SVTYPE") {
         svtype=tagelems[2];
      } else if (tagelems[1] == "DPNORM") {
         dpnorm=tagelems[2];
      } else if (tagelems[1] == "SNAME") {
         sampleid=tagelems[2];
      };
   };
   print cnvid, sampleid, svtype, dpnorm, qcnvlen, $13;
}
