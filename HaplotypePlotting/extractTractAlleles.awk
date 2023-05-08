#!/bin/awk -f
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(chrom) == 0) {
      print "Missing required variable chrom, cannot proceed." > "/dev/stderr";
      exit 2;
   };
   if (length(segcol) == 0) {
      segcol="SEGMENT";
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!(segcol in cols)) {
      print "No column named "segcol" in the input, cannot proceed." > "/dev/stderr";
      exit 3;
   };
}
NR>1{
   if ($cols["CHROM"] == chrom && $cols[segcol] == seg) {
      split($cols["ALT"], alleles, ",");
      alleles[0]=$cols["REF"];
      print $cols["CHROM"], $cols["POS"], $cols["REF"], alleles[$cols["ALLELE"]+0];
   };
}
