#!/bin/awk -f
#This script simply concatenates the core haplotypes TSVs from each
# chromosome and tract origin together, labeling them with extra
# columns indicating the S' target population and tract origin.
#This way, a downstream script can filter the .score file or the
# archaic matches TSV to only include these variants and relabel
# the tract IDs in one fell swoop.
#Required arguments:
#  pop: The name of the S' target population
#  origin: The S' tract origin
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(pop) == 0) {
      print "Missing pop variable, please set it. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(origin) == 0) {
      print "Missing origin variable, please set it. Quitting." > "/dev/stderr";
      exit 3;
   };
   firstheader=1;
}
#Only output the header line once:
/^CHROM/{
   if (firstheader) {
      print $0, "Population", "TractOrigin";
      firstheader=0;
   };
}
#Add the S' target population and tract origin columns:
!/^CHROM/{
   $3=origin"."$3;
   print $0, pop, origin;
}
