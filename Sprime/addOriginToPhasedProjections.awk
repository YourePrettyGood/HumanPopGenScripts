#!/bin/awk -f
#This script takes a match rates TSV for all S' tracts as well
# as an input BED of phased S' tract projections and reformats
# the BED into a format compatible with separateTractsByOrigin.awk.
#The phased S' tract projection BED follows the spirit of the
# BED format specification by storing custom values as tag-value
# pairs in the "name" (4th) column, so we need to convert this
# to a headered BED-like format with custom columns past the 3rd.
#These custom columns must include "Neandertal" and "Denisovan",
# which are the respective match rates as required by
# separateTractsByOrigin.awk.
#By default, we also filter the input projections by both
# ascertainable archaic sites and minimum match rates, but this
# filtering can be optionally turned off.
#Optional arguments:
#  allout:   Output all phased projections, do not filter
#  criteria: Criteria to use for filtering, currently either Browning or PFR
#            (default: PFR)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(allout) == 0) {
      allout=0;
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", allout);
   sub(/^[Nn][Oo]?$/, "0", allout);
   if (length(criteria) == 0) {
      criteria="PFR";
   };
   if (criteria != "Browning" && criteria != "PFR") {
      print "criteria "criteria" is not currently supported, please use Browning or PFR. Quitting.\n" > "/dev/stderr";
      exit 2;
   };
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#First file is the match rates file
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
#Store the ascertainable sites and match rates for each tract:
filenum==1&&FNR>1{
   #Filters to apply:
   keep=0;
   if (criteria == "Browning") {
      keep=$cols["Agood"]>=30 && $cols["Dgood"]>=30 && ($cols["Amatchrate"]>=0.3 || $cols["Dmatchrate"]>=0.3);
   } else if (criteria == "PFR") {
      keep=$cols["Ngood"]>=30 && $cols["Dgood"]>=30 && ($cols["Nmatchrate"]>=0.3 || $cols["Dmatchrate"]>=0.3);
   };
   if (allout > 0 || keep) {
      matchrates[$cols["TractID"]]=$cols["Nmatchrate"] OFS $cols["Dmatchrate"] OFS $cols["Amatchrate"] OFS $cols["Vmatchrate"] OFS $cols["Cmatchrate"];
   };
}
#Second file is the phased projections BED
#Optionally, any further files are also processed as phased
# projections BEDs with the same tag-value pairs as the
# second file.
#Compose the header, determining extra columns from the tags
# of the first line of the input:
filenum==2&&FNR==1{
   outcolnames[1]="#Chromosome";
   outcolnames[2]="Start";
   outcolnames[3]="End";
   n_tags=split($4, tags, ";");
   for (i=1; i<=n_tags; i++) {
      split(tags[i], tagelems, "=");
      outcolnames[3+i]=tagelems[1];
   };
   outcolnames[3+n_tags+1]="Neandertal";
   outcolnames[3+n_tags+2]="Denisovan";
   outcolnames[3+n_tags+3]="AltaiNea";
   outcolnames[3+n_tags+4]="VindijaNea";
   outcolnames[3+n_tags+5]="ChagyrskayaNea";
   n_outcols=3+n_tags+5;
   printf "%s", outcolnames[1];
   for (i=2; i<=n_outcols; i++) {
      printf "%s%s", OFS, outcolnames[i];
   };
   printf "\n";
}
#Process all projections after adding the header:
filenum>=2{
   delete linevals;
   linevals["#Chromosome"]=$1;
   linevals["Start"]=$2;
   linevals["End"]=$3;
   n_tags=split($4, tags, ";");
   for (i=1; i<=n_tags; i++) {
      split(tags[i], tagelems, "=");
      linevals[tagelems[1]]=tagelems[2];
   };
   if (linevals["TractID"] in matchrates) {
      printf "%s", linevals[outcolnames[1]];
      for (i=2; i<=n_tags+3; i++) {
         printf "%s%s", OFS, linevals[outcolnames[i]];
      };
      printf "%s%s", OFS, matchrates[linevals["TractID"]];
      printf "\n";
   };
}
