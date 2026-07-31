#!/bin/awk -f
#This takes the output of bcftools query -H -f '%ID\t%REF\t%ALT\t%INFO/AF_...'
# where '%INFO/AF_...' is an arbitrarily long string specifying other tags
# such as allele frequencies you want to transfer to the output of this
# script.
#The general idea is that we can keep track of allele frequencies for SVs
# from the VCF input to bcftools query as well as sequence information
# necessary for downstream analyses like CADD-SV.
#Assuming the VCF came from PanGenie, the ID will contain CHROM, POS, SVTYPE,
# and SVLEN, as well as a sequence of pangenome graph nodes and edges defining
# the variant. We parse these apart in order to filter for SVs only and output
# SVTYPE and SVLEN.
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
   if (!($1 ~ /^#/)) {
      print "Input file appears to be missing a header. Make sure the bcftools query call includes -H to print a header. Quitting." > "/dev/stderr";
      exit 2;
   };
   if ($1 != "#[1]ID") {
      print "First column of input does not appear to be ID. Quitting." > "/dev/stderr";
      exit 3;
   };
   if ($2 != "[2]REF") {
      print "Second column of input does not appear to be REF. Quitting." > "/dev/stderr";
      exit 4;
   };
   if ($3 != "[3]ALT") {
      print "Third column of input does not appear to be ALT. Quitting." > "/dev/stderr";
      exit 5;
   };
   printf "CHROM\tPOS\tSVTYPE\tSVLEN";
   for (i=2; i<=NF; i++) {
      sub("[[0-9]+]", "", $i);
      printf "\t%s", $i;
   };
   printf "\n";
}
NR>1{
   n=split($1, a, "-");
   chrom=a[1];
   pos=a[2];
   svtype=a[3];
   svlen=a[5];
   if (svlen >= 50) {
      printf "%s\t%s\t%s\t%s", chrom, pos, svtype, svlen;
      for (i=2; i<=NF; i++) {
         printf "\t%s", $i;
      };
      printf "\n";
   };
}
