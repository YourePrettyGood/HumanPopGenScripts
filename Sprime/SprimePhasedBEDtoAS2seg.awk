#!/bin/awk -f
#This script takes in the tract match rates TSV and the
# BED produced by SprimeTractBED.awk with phased=1 in order
# to generate an approximation of the .seg file produced
# by ArchaicSeeker2.0.
#The idea here is to generate a file compatible with
# getAS2Seg so that we can run MultiWaver2.1 on Sprime
# tracts.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(criteria) == 0) {
      criteria="PFR";
   } else if (criteria != "Browning" && criteria != "PFR") {
      print "Invalid match criteria name ("criteria"), please specify either Browning or PFR" > "/dev/stderr";
      exit 2;
   };
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#First input file is the decompressed Sprime match rates TSV.
#We primarily just need the following columns to establish
# tract origin:
# - TractID
# - Ngood
# - Agood
# - Dgood
# - Amatchrate
# - Nmatchrate
# - Dmatchrate
#Construct a map from column names to column indices:
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      mrcols[$i]=i;
   };
}
#Classify tract origins for the last column of the .seg output:
#The Browning thresholds here match the paper for NEA and DEN,
# but Ambiguous isn't ever established in the paper, so we set
# a rough threshold here.
filenum==1&&FNR>1{
   Mthresh["Browning"]=$mrcols["Agood"] >= 30 && $mrcols["Dgood"] >= 30;
   Nthresh["Browning"]=$mrcols["Amatchrate"] >= 0.6 && $mrcols["Dmatchrate"] <= 0.4;
   Dthresh["Browning"]=$mrcols["Amatchrate"] <= 0.3 && $mrcols["Dmatchrate"] >= 0.4;
   Athresh["Browning"]=$mrcols["Amatchrate"] > 0.3 && $mrcols["Dmatchrate"] > 0.4;
   Mthresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30;
   Nthresh["PFR"]=$mrcols["Nmatchrate"] >= 0.3 && $mrcols["Dmatchrate"] <= 0.3;
   Dthresh["PFR"]=$mrcols["Nmatchrate"] <= 0.3 && $mrcols["Dmatchrate"] >= 0.3;
   Athresh["PFR"]=$mrcols["Nmatchrate"] > 0.3 && $mrcols["Dmatchrate"] > 0.3;
   if (Mthresh[criteria]) {
      if (Nthresh[criteria]) {
         origin[$mrcols["TractID"]]="Neandertal";
      } else if (Dthresh[criteria]) {
         origin[$mrcols["TractID"]]="Denisovan";
      } else if (Athresh[criteria]) {
         origin[$mrcols["TractID"]]="Ambiguous";
      } else {
         origin[$mrcols["TractID"]]="BadMatch";
      };
   } else {
      origin[$mrcols["TractID"]]="TooMuchMissing";
   };
}
#The second input file is the output of SprimeTractBED.awk with phased=1.
#First output the header line for the .seg:
#.seg format: Haplotype ID, CHROM, START, END, MAPSTART, MAPEND, ORIGIN[, ...]
filenum==2&&FNR==1{
   print "ID", "Contig", "Start(bp)", "End(bp)", "Start(cM)", "End(cM)", "BestMatchedPop";
}
#We extract the tract ID and haplotype ID from the fourth column and
# output a .seg line with start adjusted to be 1-based, as well as the
# tract origin we calculated above:
filenum==2{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], tagelems, "=");
      if (tagelems[1] == "TractID") {
         tractid=tagelems[2];
      } else if (tagelems[1] == "Haplotype") {
         id=tagelems[2];
      };
   };
   print id, $1, $2+1, $3, "-", "-", origin[tractid];
}
