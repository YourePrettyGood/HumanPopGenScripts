#!/bin/awk -f
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(source) == 0) {
      print "Please specify an archaic source to select tracts from." >> "/dev/stderr";
      exit 2;
   };
   if (length(only_matches) == 0) {
      only_matches=0;
   };
   if (length(criteria) == 0) {
      criteria="Browning";
   } else if (criteria != "Browning" && criteria != "PFR") {
      print "Invalid match criteria name ("criteria"), please specify either Browning or PFR" > "/dev/stderr";
      exit 3;
   };
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#First file is the Sprime match rates file we use to select tracts
# from a particular archaic source:
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      mrcols[$i]=i;
   };
}
#Retain tracts that come from whichever source is selected:
filenum==1&&FNR>1{
   Nthresh["Browning"]=$mrcols["Agood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Amatchrate"] >= 0.6 && $mrcols["Dmatchrate"] <= 0.4;
   Dthresh["Browning"]=$mrcols["Agood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Amatchrate"] <= 0.3 && $mrcols["Dmatchrate"] >= 0.4;
   Nthresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Nmatchrate"] >= 0.3 && $mrcols["Dmatchrate"] <= 0.3;
   Dthresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Nmatchrate"] <= 0.3 && $mrcols["Dmatchrate"] >= 0.3;
   if ((source == "Neandertal" && Nthresh[criteria]) || (source == "Denisovan" && Dthresh[criteria])) {
      keep[$mrcols["TractID"]]=1;
   };
}
#Second file is the Sprime matches file (precursor to the match rates
# file) so we can extract Sprime sites from the selected tracts:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      mcols[$i]=i;
   };
}
#Print out sites if they're in the selected tracts, optionally
# filtering sites by whether or not they match the source
# archaic:
filenum==2&&FNR>1{
   if ($mcols["TractID"] in keep) {
      if (only_matches && $mcols[source"Match"] == "match") {
         if (length(debug) > 0 && debug >= 2) {
            print $0;
         } else {
            print $mcols["CHROM"], $mcols["POS"];
         };
      } else if (!only_matches) {
         if (length(debug) > 0 && debug >= 2) {
            print $0;
         } else {
            print $mcols["CHROM"], $mcols["POS"];
         };
      };
   };
}
