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
      } else if (tagelems[1] == "RefLength") {
         rlen=tagelems[2];
      };
   };
   if (!approxlen && prevqname == qname && prevchr == $1 && $2 < prevend) {
      #This case covers contained or overlapping alignments, respectively:
      if ($3-prevend <= 0) {
         #If the current alignment is contained in the previous one, the
         # additional length contributed is 0, though the end of the
         # current alignment minus that of the previous alignment is
         # zero or negative.
         addllen=0;
      } else {
         #If the current alignment extends beyond the previous one, add
         # the extension.
         addllen=$3-prevend;
      };
   } else {
      #This case covers either if we want approximate length, or if
      # alignments are abutting or disjoint:
      addllen=$3-$2;
   };
   if (addllen < 0) {
      print $0, prevchr, prevstart, prevend, addllen, prevqname > "/dev/stderr";
   };
   alnlensum[qname,$1]+=addllen;
   if (unweightedorient) {
      #Just count up alignments as votes:
      alnorient[qname,$1,$6]++;
   } else {
      #Weigh each alignment's vote by the additional length it contributes:
      alnorient[qname,$1,$6]+=addllen;
   };
   qlens[qname]=qlen;
   rlens[$1]=rlen;
   prevchr=$1;
   prevend=$3;
   prevstart=$2;
   prevqname=qname;
}
END{
   PROCINFO["sorted_in"]="@val_num_desc";
   #Old notes prior to introducing approxlen and unweightedorient:
   #Note that the alignment length calculated here is based on the
   # reference sequence portion of the alignment, so it's possible
   # to have > 100% query coverage using this calculation.
   #Furthermore, I don't adjust for overlapping alignments, so that
   # can also produce query coverage > 100%.
   #Orientation is determined using majority vote, not weighting
   # by alignment length (which would probably be better).
   if (approxlen) {
      print "Query", "Reference", "~QueryCov", "~RefCov", "Strand", "~AlnRefLen", "QueryLen", "RefLen", "+Alns", "-Alns";
   } else if (unweightedorient) {
      print "Query", "Reference", "QueryCov", "RefCov", "Strand", "AlnRefLen", "QueryLen", "RefLen", "+Alns", "-Alns";
   } else {
      print "Query", "Reference", "QueryCov", "RefCov", "Strand", "AlnRefLen", "QueryLen", "RefLen", "+AlnLen", "-AlnLen";
   };
   for (qrmatch in alnlensum) {
      split(qrmatch, queryref, SUBSEP);
      qname=queryref[1];
      rname=queryref[2];
      asum=alnlensum[qrmatch];
      qlen=qlens[qname];
      rlen=rlens[rname];
      qcov=asum*100/qlen;
      rcov=asum*100/rlen;
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
      print qname, rname, qcov, rcov, orient, asum, qlen, rlen, posorient, negorient;
   };
}
