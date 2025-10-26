#!/bin/awk -f
#This script takes in a sample metadata file produced by ts_to_metadata.awk
# and generates a TSV suitable for use as the -G option of bcftools +split
# in order to split a VCF into one VCF containing only modern samples and
# the other VCF containing only archaic samples.
#Samples are sorted into these two groups based on the value of the second
# column of the metadata file (assuming sample ID is in the first column
# and that there is a header line) matching one of the population IDs listed
# in either the "arc" or "modern" argument lists.
#Required arguments:
# prefix: Prefix for the VCFs to be output by bcftools +split
# popcol: Column name for population IDs
#         (default: Region)
# arc:    Comma-separated list of population IDs of archaic samples
#         in the sample metadata file provided
# modern: Comma-separated list of population IDs of modern samples
#         in the sample metadata file provided
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(prefix) == 0) {
      print "No prefix for output VCFs provided. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(popcol) == 0) {
      popcol="Region";
   };
   if (length(arc) == 0) {
      print "No arc provided, please provide a comma-separated list of population names that should be mapped to the archaic output VCF. Quitting." > "/dev/stderr";
      exit 3;
   };
   n_arc=split(arc, arcarr, ",");
   for (i=1; i<=n_arc; i++) {
      archash[arcarr[i]]=i;
   };
   if (length(modern) == 0) {
      print "No modern provided, please provide a comma-separated list of population names that should be mapped to the modern output VCF. Quitting." > "/dev/stderr";
      exit 4;
   };
   n_modern=split(modern, modernarr, ",");
   for (i=1; i<=n_modern; i++) {
      modernhash[modernarr[i]]=i;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!(popcol in cols)) {
      print "popcol ("popcol") not found in metadata file header, cannot proceed." > "/dev/stderr";
      exit 5;
   };
}
NR>1{
   if ($cols[popcol] in archash) {
      print $cols["SampleID"], "-", prefix"_archaic.vcf.gz";
   } else if ($cols[popcol] in modernhash) {
      print $cols["SampleID"], "-", prefix"_modern.vcf.gz";
   };
}
