#!/bin/awk -f
#This script takes the concatenated core haplotypes as
# well as the S' .score file or a file with those columns
# and more (e.g. the archaic matches TSV), filters the
# sites to only retain those found in the core haplotypes,
# adjusts their SEGMENT ID to the cluster ID of the core
# haplotype, and adds extra columns for metadata and
# backtracing.
#The output of this script should be directly compatible
# with SprimeArchaicAF.awk and SprimeArchaicMatchRate.awk.
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#First file is the concatenated core haplotypes file:
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      chcols[$i]=i;
   };
}
#Store the cluster ID as the new SEGMENT, as well as the
# S' target population and the tract origin:
filenum==1&&FNR>1{
   tract[$chcols["CHROM"],$chcols["POS"]]=$chcols["ClusterID"];
   spop[$chcols["CHROM"],$chcols["POS"]]=$chcols["Population"];
   orig[$chcols["CHROM"],$chcols["POS"]]=$chcols["TractOrigin"];
}
#Second file is the S' .score file or the archaic matches
# TSV or anything that retains the S' .score file's columns:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      scols[$i]=i;
   };
   #Print the modified header:
   print $0, "TargetPopulation", "TractOrigin", "OriginalSEGMENT";
}
#Only output sites that are found in the core haplotypes,
# and make sure to relabel their SEGMENT based on the
# cluster ID of the core haplotype:
filenum==2&&FNR>1{
   site=$scols["CHROM"] SUBSEP $scols["POS"];
   if (site in tract) {
      origseg=$scols["SEGMENT"];
      $scols["SEGMENT"]=tract[site];
      print $0, spop[site], orig[site], origseg;
   };
}
