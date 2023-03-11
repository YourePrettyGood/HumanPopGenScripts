#!/bin/awk -f
#This script takes a modified .score or matches file produced by
# applyCoreHaplotypesToScore.awk and identifies the correspondence
# between Sprime TractIDs (taken as
# [TargetPopulation]_[CHROM]_[OriginalSEGMENT]) and the core
# haplotype ID (taken as [TargetPopulation]_[CHROM]_[SEGMENT]),
# deduplicating along the way (since the input file contains one
# line per Sprime variant in each core haplotype).
#Options:
# header: (optional) A flag indicating whether or not to print the header
#                    line (default: don't print header)
#The "header" flag facilitates combining the results across chromosomes and
# populations.
BEGIN{
   FS="\t";
   OFS=FS;
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", header);
   sub(/^[Nn][Oo]?$/, "0", header);
   if (length(header) > 0 && header > 0) {
      print "CorehapID", "TractID";
   };
}
#Parse the header line for column names:
NR==1{
   for (i=1; i<=NF; i++) {
     cols[$i]=i;
   };
}
#Now construct the core haplotype ID for each line and output
# unique mappings:
NR>1{
   tractid=$cols["TargetPopulation"]"_"$cols["CHROM"]"_"$cols["OriginalSEGMENT"];
   corehapid=$cols["TargetPopulation"]"_"$cols["CHROM"]"_"$cols["SEGMENT"];
   if (!((corehapid,tractid) in dedup)) {
      print corehapid, tractid;
      dedup[corehapid,tractid]++;
   };
}
