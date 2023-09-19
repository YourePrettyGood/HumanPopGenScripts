#!/bin/awk -f
#This script filters a BED of Sprime tract projections to
# only retain projections for tracts with TractID tag values
# found in the input list.
#The first input file should be a tab-separated file with
# a column named "TractID" containing the Sprime tract IDs
# to be retained.
#The second input file is the BED of Sprime tract
# projections to be filtered.
BEGIN{
   FS="\t";
   OFS=FS;
}
FNR==NR&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
FNR==NR&&FNR>1{
   keep[$cols["TractID"]]++;
}
FNR<NR{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], tagelems, "=");
      if (tagelems[1] == "TractID") {
         if (tagelems[2] in keep) {
            print;
         };
      };
   };
}
