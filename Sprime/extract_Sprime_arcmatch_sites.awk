#!/bin/awk -f
#This script identifies a set of variant sites from S' tracts of
# a particular origin (as determined by match rate criteria)
# for downstream use in determining core haplotypes of S' tracts.
#That sounds a bit vague/abstract, but the general idea is that
# this is a generalization of the first filter applied in
# Gittelman et al. 2016 for identification of core haplotypes
# from S* tracts. In Gittelman et al., S* tag SNPs from
# windows with adequate density are filtered by whether or not
# they match the Altai Neanderthal. One of S*'s steps already
# classifies each tract based on it's archaic origin, so we
# adapt both of these filters to S' variants.
#Thus, we take in the matches and match rates files output by
# post-processing S' .score files, classify each S' tract
# based on a set of match rate criteria (currently limited to
# "Browning" and "PFR"), and only output the chromosome and
# position of S' variants in these tracts, optionally further
# filtering by whether or not they match a particular archaic
# source. For Neandertal tracts, we check if the variant matches
# any of the 3 high-depth Neanderthals; for Denisovan tracts,
# we check if the variant matches the Altai Denisovan; and for
# Ambiguous tracts we check if the variant matches any of the
# 3 high-depth Neanderthals or the Altai Denisovan.
#The resulting TSV of variant sites can then be fed into the
# -R and -T arguments of `bcftools view` to then feed into
# vcftools for calculation of pairwise r^2 between all of these
# variants, filtering for r^2 >= 0.3 as per Gittelman et al. 2016
# and Vernot et al. 2016.
#Required arguments:
#  source: Either Neandertal, Denisovan, or Ambiguous
#Optional arguments:
#  criteria: Which match rate criteria to use (either Browning or PFR)
#            (default: Browning)
#  only_matches: Whether to only output sites that match the appropriate
#                archaic, or output all qualifying sites
#                (default: output all qualifying sites)
#First input file should be the decompressed match rates file
# produced by SprimeArchaicMatchRate.awk
#Second input file should be the decompressed matches file
# produced by archaicMatchSprime.awk
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
   Athresh["Browning"]=0;
   Nthresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Nmatchrate"] >= 0.3 && $mrcols["Dmatchrate"] <= 0.3;
   Dthresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Nmatchrate"] <= 0.3 && $mrcols["Dmatchrate"] >= 0.3;
   Athresh["PFR"]=$mrcols["Ngood"] >= 30 && $mrcols["Dgood"] >= 30 && $mrcols["Nmatchrate"] >= 0.3 && $mrcols["Dmatchrate"] >= 0.3;
   if ((source == "Neandertal" && Nthresh[criteria]) || (source == "Denisovan" && Dthresh[criteria]) || (source == "Ambiguous" && Athresh[criteria])) {
      keep[$mrcols["TractID"]]=1;
   };
}
#Second file is the Sprime matches file (precursor to the match rates
# file) so we can extract Sprime sites from the selected tracts:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      mcols[$i]=i;
   };
   #Also handle the ambiguous tract origin case:
   if (source == "Ambiguous") {
      source="Neandertal";
      altsource="Denisovan";
   } else {
      altsource=source;
   };
}
#Print out sites if they're in the selected tracts, optionally
# filtering sites by whether or not they match the source
# archaic:
filenum==2&&FNR>1{
   if ($mcols["TractID"] in keep) {
      if (only_matches && ($mcols[source"Match"] == "match" || $mcols[altsource"Match"] == "match")) {
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
