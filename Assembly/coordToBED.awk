#!/bin/awk -f
#This script does a simple conversion from the output of MUMmer's
# show-coords tool with options -cdHloqT to a BED6 file.
#The BED6 is from the point of view of the reference by default,
# as this enables downstream selection of the best syntenic hit
# for each query contig/scaffold. The score column (#5) contains
# the percent identity of the alignment, and the name column (#4)
# contains a GFF3-like tag string with the query contig/scaffold
# name and the overall length of that query contig.
#2024/09/06: Added an option to process PAF from minimap2/minigraph
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(format)==0) {
      format="coord";
   } else if (format ~ /^[Mm][Uu][Mm][Mm][Ee][Rr]|^[Cc][Oo][Oo][Rr][Dd]$/) {
      format="coord";
   } else if (format ~ /^[Mm][Ii][Nn][Ii]|^[Pp][Aa][Ff]$/) {
      format="paf";
   } else {
      print "Unrecognized format "format > "/dev/stderr";
      exit 2;
   };
}
{
   if (format == "coord") {
      rname=$14;
      qname=$15;
      pident=$7;
      strand=$13=="-1"?"-":"+";
      rlen=$8;
      qlen=$9;
      #Coordinates here should be 1-based:
      qstart=$3;
      qend=$4;
      rstart=$1;
      rend=$2;
   } else if (format == "paf") {
      rname=$6;
      #Make sure to skip unmapped contigs:
      if (rname == "*") {
         next;
      };
      qname=$1;
      if ($11 > 0) {
         pident=$10*100/$11;
      } else {
         pident="NA";
      };
      strand=$5;
      rlen=$7;
      qlen=$2;
      #Coordinates here should be 1-based, so add 1:
      qstart=$3+1;
      qend=$4+1;
      rstart=$8+1;
      rend=$9+1;
   } else {
      print "Unknown input format "format > "/dev/stderr";
      exit 2;
   };
   if (length(flip)>0) {
      if (strand == "-") {
         astart=qend-1;
         aend=qstart;
      } else {
         astart=qstart-1;
         aend=qend;
      };
      print qname, astart, aend, "QueryName="rname";QueryLength="rlen";RefLength="qlen, pident, strand;
   } else {
      astart=rstart-1;
      aend=rend;
      print rname, astart, aend, "QueryName="qname";QueryLength="qlen";RefLength="rlen, pident, strand;
   };
}
