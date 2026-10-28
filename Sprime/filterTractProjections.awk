#!/bin/awk -f
#This script filters a BED of Sprime tract projections to
# only retain projections for tracts with TractID tag values
# found in the input list.
#The first input file should be a tab-separated file with
# a column named "TractID" containing the Sprime tract IDs
# to be retained.
#The second input file is the BED of Sprime tract
# projections to be filtered.
#2026/10/14: Generalized to allow for filtering by any
# user-specified tag, not just TractID.
# TractID is now the default for backwards compatibility.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(tagname) == 0) {
      tagname="TractID";
   };
}
FNR==NR&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!(tagname in cols)) {
      print "Did not find "tagname" tag in the first input file. Quitting." > "/dev/stderr";
      exit 2;
   };
}
FNR==NR&&FNR>1{
   keep[$cols[tagname]]++;
}
FNR<NR{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], tagelems, "=");
      if (tagelems[1] == tagname) {
         if (tagelems[2] in keep) {
            print;
         };
      };
   };
}
