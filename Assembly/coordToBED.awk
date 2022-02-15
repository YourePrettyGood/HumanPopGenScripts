#!/bin/awk -f
#This script does a simple conversion from the output of MUMmer's
# show-coords tool with options -cdHloqT to a BED6 file.
#The BED6 is from the point of view of the reference by default,
# as this enables downstream selection of the best syntenic hit
# for each query contig/scaffold. The score column (#5) contains
# the percent identity of the alignment, and the name column (#4)
# contains a GFF3-like tag string with the query contig/scaffold
# name and the overall length of that query contig.
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   rname=$14;
   qname=$15;
   pident=$7;
   strand=$13=="-1"?"-":"+";
   rlen=$8;
   qlen=$9;
   if (length(flip)>0) {
      if (strand == "-") {
         astart=$4-1;
         aend=$3;
      } else {
         astart=$3-1;
         aend=$4;
      };
      print qname, astart, aend, "QueryName="rname";QueryLength="rlen";RefLength="qlen, pident, strand;
   } else {
      astart=$1-1;
      aend=$2;
      print rname, astart, aend, "QueryName="qname";QueryLength="qlen";RefLength="rlen, pident, strand;
   };
}
