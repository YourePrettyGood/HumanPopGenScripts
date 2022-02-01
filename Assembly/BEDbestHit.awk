#!/bin/awk -f
#This script takes the unmerged BED of nucmer alignments and identifies
# the best hit for each query contig/scaffold based on approximate
# summed query coverage, and picks the best orientation based on
# majority vote across the component alignments.
#A better estimator of orientation probably would weight the votes by
# the length of each alignment.
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], tagelems, "=");
      if (tagelems[1] == "QueryName") {
         qname=tagelems[2];
      } else if (tagelems[1] == "QueryLength") {
         qlen=tagelems[2];
      };
   };
   alnlensum[qname,$1]+=$3-$2;
   alnorient[qname,$1,$6]++;
   qlens[qname]=qlen;
}
END{
   PROCINFO["sorted_in"]="@val_num_desc";
   #Note that the alignment length calculated here is based on the
   # reference sequence portion of the alignment, so it's possible
   # to have > 100% query coverage using this calculation.
   #Furthermore, I don't adjust for overlapping alignments, so that
   # can also produce query coverage > 100%.
   #Orientation is determined using majority vote, not weighting
   # by alignment length (which would probably be better).
   print "Query", "Reference", "~QueryCov", "~AlnRefLen", "QueryLen", "MajorityStrand", "+Alns", "-Alns";
   for (qrmatch in alnlensum) {
      split(qrmatch, queryref, SUBSEP);
      qname=queryref[1];
      rname=queryref[2];
      asum=alnlensum[qrmatch];
      qlen=qlens[qname];
      qcov=asum*100/qlen;
      posorient=length(alnorient[qname,rname,"+"])>0?alnorient[qname,rname,"+"]:0;
      negorient=length(alnorient[qname,rname,"-"])>0?alnorient[qname,rname,"-"]:0;
      if (posorient + negorient > 0) {
         if (posorient > negorient) {
            orient="+";
         } else if (negorient > posorient) {
            orient="-";
         } else {
            orient="?";
         };
      } else {
         orient="?";
      };
      print qname, rname, qcov, asum, qlen, orient, posorient, negorient;
#      print queryref[1], queryref[2], alnlensum[qrmatch]*100/qlens[queryref[1]], alnlensum[qrmatch], qlens[queryref[1]], alnorient[qrmatch,"-"]>alnorient[qrmatch,"+"]?"-":"+", alnorient[qrmatch,"+"]>0?alnorient[qrmatch,"+"]:0, alnorient[qrmatch,"-"]>0?alnorient[qrmatch,"-"]:0;
   };
}
