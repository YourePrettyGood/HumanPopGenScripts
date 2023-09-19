#!/bin/awk -f
#This script splits an input BED of S' tracts or core haplotypes
# into separate files based on tract origin (as determined by
# match rate criteria).
#Expected columns include "Neandertal" and "Denisovan",
# which are the match rates to all 3 Neandertals and the
# Altai Denisovan, respectively.
#Output files will have names of the format:
# [prefix]_[origin]_[suffix]
# where [origin] can be Neandertal, Denisovan, or Ambiguous
#Required arguments:
#  prefix:   Prefix for the output file names
#  suffix:   Suffix for the output file names
#Optional arguments:
#  criteria: Which match rate criteria to use (either Browning or PFR)
#            (default: PFR)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(prefix) == 0) {
      print "Please specify an output file prefix." >> "/dev/stderr";
      exit 2;
   };
   if (length(suffix) == 0) {
      print "Please specify an output file suffix." >> "/dev/stderr";
      exit 3;
   };
   if (length(criteria) == 0) {
      criteria="PFR";
   } else if (criteria != "Browning" && criteria != "PFR") {
      print "Invalid match criteria name ("criteria"), please specify either Browning or PFR" > "/dev/stderr";
      exit 4;
   };
   origin[1]="Ambiguous";
   origin[2]="Denisovan";
   origin[3]="Neandertal";
   PROCINFO["sorted_in"]="@ind_num_asc";
}
#Keep track of the column names and output the header:
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   for (i in origin) {
      print > prefix"_"origin[i]"_"suffix;
   };
}
#For each line, evaluate the match rate criteria and classify the tract:
NR>1{
   outfn="";
   thresh["Neandertal","Browning"]=$cols["AltaiNea"] >= 0.6 && $cols["Denisovan"] <= 0.4;
   thresh["Denisovan","Browning"]=$cols["AltaiNea"] <= 0.3 && $cols["Denisovan"] >= 0.4;
   thresh["Ambiguous","Browning"]=$cols["AltaiNea"] > 0.3 && $cols["Denisovan"] > 0.4;
   thresh["Neandertal","PFR"]=$cols["Neandertal"] >= 0.3 && $cols["Denisovan"] <= 0.3;
   thresh["Denisovan","PFR"]=$cols["Neandertal"] <= 0.3 && $cols["Denisovan"] >= 0.3;
   thresh["Ambiguous","PFR"]=$cols["Neandertal"] > 0.3 && $cols["Denisovan"] > 0.3;
   for (i in origin) {
      if (thresh[origin[i],criteria]) {
         outfn=prefix"_"origin[i]"_"suffix;
      };
   };
   if (outfn != "") {
      print > outfn;
   };
}
