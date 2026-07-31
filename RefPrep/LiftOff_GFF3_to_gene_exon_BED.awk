#!/bin/awk -f
#This is a very quick and rough script to convert LiftOff-based GFF3 annotations
# into BED format, specifically focusing on gene and exon features.
#It's mainly intended for the LiftOff annotations of T2T-CHM13 v2.0.
BEGIN{
   FS="\t";
   OFS=FS;
}
!/^#/{
   if ($3 == "gene") {
      #Parse the tag-value pairs for gene_name and gene_biotype:
      n=split($9, a, ";");
      for (i=1; i<=n; i++) {
         m=split(a[i], b, "=");
         if (b[1] == "gene_name") {
            id=b[2];
         } else if (b[1] == "gene_biotype") {
            type=b[2];
         };
      };
      print $1, $4-1, $5, "FeatureType=gene;ID="id";Biotype="type;
   } else if ($3 == "exon") {
      #Parse the tag-value pairs for gene (name of the gene the exon belongs to):
      n=split($9, a, ";");
      for (i=1; i<=n; i++) {
         m=split(a[i], b, "=");
         if (b[1] == "gene") {
            id=b[2];
         };
      };
      print $1, $4-1, $5, "FeatureType=exon;ID="id;
   };
}
