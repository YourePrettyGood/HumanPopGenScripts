#!/bin/awk -f
#This script takes the output of the R script that filters PanGenie-based SV selection scan results
# based on significance from simulations and reformats the results into BED-like records for each
# SV.
#Use this as the -a input for `bedtools intersect -wao` to do some very rough annotation of
# SV consequences.
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
   #Store the column indices by name to deal with swapped column order:
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!("CHROM" in cols) || !("POS" in cols) || !("SVTYPE" in cols) || !("SVLEN" in cols)) {
      print "Missing one or more of CHROM, POS, SVTYPE, and SVLEN columns from header of input. Quitting." > "/dev/stderr";
      exit 2;
   };
}
NR>1{
   #Reformat the SVs into a BED-style format, with GFF3-style tag-value pairs in the fourth column:
   #Also, we adjust the BED interval to represent deletions, but not insertions and complex SVs.
   if ($cols["SVTYPE"] == "DEL") {
      print $cols["CHROM"], $cols["POS"]-1, $cols["POS"]+$cols["SVLEN"]-1, "SVTYPE="$cols["SVTYPE"]";SVLEN="$cols["SVLEN"];
   } else {
      print $cols["CHROM"], $cols["POS"]-1, $cols["POS"], "SVTYPE="$cols["SVTYPE"]";SVLEN="$cols["SVLEN"];
   };
}
