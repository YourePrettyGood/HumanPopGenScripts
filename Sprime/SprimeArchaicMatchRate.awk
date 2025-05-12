#!/bin/awk -f
#This script takes the output of archaicMatchSprime.awk and summarizes the
# match indicators into match rates and counts of ascertainable sites.
#Optional arguments:
# debug:   Output debug info to STDERR
# species: Comma-separated list of archaic hominin groups to summarize
#          Default: Neandertal,Denisovan
# arcs:    Comma-separated list of archaic hominin samples to summarize
#          Default: Altai,Vindija,Chagyrskaya
# short:   Comma-separated list of single characters representing
#          groups,arcs
#          Default: first letter of each
BEGIN{
   FS="\t";
   OFS=FS;
   #Extract the archaic hominin groups and individuals:
   if (length(species) == 0) {
      species="Neandertal,Denisovan";
   };
   n_spp=split(species, spp, ",");
   if (length(arcs) == 0) {
      arcs="Altai,Vindija,Chagyrskaya";
   };
   n_arcs=split(arcs, arcindivs, ",");
   #Combine the two and set the abbreviations/shortgroups:
   for (i=1; i<=n_spp; i++) {
      groups[i]=spp[i];
      if (length(short) == 0) {
         shortgroups[i]=substr(spp[i], 1, 1);
      };
      n_groups+=1;
   };
   for (i=1; i<=n_arcs; i++) {
      groups[n_spp+i]=arcindivs[i];
      if (length(short) == 0) {
         shortgroups[n_spp+i]=substr(arcindivs[i], 1, 1);
      };
      n_groups+=1;
   };
   if (length(short) > 0) {
      n_short=split(short, shortgroups, ",");
      if (n_short != n_groups) {
         print "Number of group/individual shortnames/abbreviations in shorts does not match the number of groups+individuals. Quitting." > "/dev/stderr";
         exit 16;
      };
   };
#   #Archaic groups are hard-coded for now:
#   groups[1]="Neandertal";
#   shortgroups[1]="N";
#   groups[2]="Denisovan";
#   shortgroups[2]="D";
#   groups[3]="Altai";
#   shortgroups[3]="A";
#   groups[4]="Vindija";
#   shortgroups[4]="V";
#   groups[5]="Chagyrskaya"
#   shortgroups[5]="C";
#   ngroups=5;
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   problem=0;
   if (!("CHROM" in cols)) {
      print "CHROM column missing from input." > "/dev/stderr";
      problem+=1;
   };
   if (!("POS" in cols)) {
      print "POS column missing from input." > "/dev/stderr";
      problem+=2;
   };
   if (!("TractID" in cols)) {
      print "TractID column missing from input." > "/dev/stderr";
      problem+=4;
   };
   if (!("SCORE" in cols)) {
      print "SCORE column missing from input." > "/dev/stderr";
      problem+=8;
   };
   if (problem > 0) {
      exit problem;
   };
}
NR>1{
   #For each tract, keep track of the interval (start=min(POS), end=max(POS)):
   #Also keep track of the score for diagnostics/filtering.
   if ($cols["TractID"] in tractchrom) {
      if ($cols["POS"] < tractstart[$cols["TractID"]]) {
         tractstart[$cols["TractID"]]=$cols["POS"];
      };
      if ($cols["POS"] > tractend[$cols["TractID"]]) {
         tractend[$cols["TractID"]]=$cols["POS"];
      };
      if ($cols["SCORE"] != tractscore[$cols["TractID"]]) {
         print "Scores for "$cols["TractID"]" are not all identical" > "/dev/stderr";
         print $cols["SCORE"]" != "tractscore[$cols["TractID"]] > "/dev/stderr";
      };
   } else {
      tractchrom[$cols["TractID"]]=$cols["CHROM"];
      tractstart[$cols["TractID"]]=$cols["POS"];
      tractend[$cols["TractID"]]=$cols["POS"];
      tractscore[$cols["TractID"]]=$cols["SCORE"];
   };
   #Also keep track of the number of sites in the tract:
   #Note: We don't deduplicate this count, although in principle there shouldn't
   # be duplicate sites.
   tractsnplen[$cols["TractID"]]++;
   #Keep track of match counts and ascertainable counts:
   for (g=1; g<=n_groups; g++) {
      if ($cols[groups[g]"Match"] == "match") {
         tractmatches[$cols["TractID"],groups[g]]+=1;
         tractascertainable[$cols["TractID"],groups[g]]+=1;
      } else if ($cols[groups[g]"Match"] == "mismatch") {
         tractascertainable[$cols["TractID"],groups[g]]+=1;
      };
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   printf "CHROM\tSTART\tEND\tTractID\tSNPLEN\tSCORE";
   for (g=1; g<=n_groups; g++) {
      printf "\t%smatchrate\t%sgood", shortgroups[g], shortgroups[g];
   };
   printf "\n";
   for (t in tractchrom) {
      printf "%s\t%s\t%s\t%s\t%i\t%s", tractchrom[t], tractstart[t], tractend[t], t, tractsnplen[t], tractscore[t];
      for (g=1; g<=n_groups; g++) {
         if (tractascertainable[t,groups[g]] > 0) {
            printf "\t%f\t%i", tractmatches[t,groups[g]]/tractascertainable[t,groups[g]], tractascertainable[t,groups[g]];
         } else {
            printf "\tNA\t0";
         };
      };
      printf "\n";
   };
}
