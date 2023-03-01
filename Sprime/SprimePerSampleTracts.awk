#!/bin/awk -f
#This script takes the S' .score file and the output of bcftools query -H
# -f '%CHROM:%POS[\t%GT]\n' and outputs sites in each sample where the
# sample has at least one archaic allele (indicating the state of the
# sample as heterozygous or homozygous).
#Required:
# spop:              The ID of the Sprime target population, used as the
#                    prefix for the tract ID
#                    (default for backwards compatibility: empty)
#Options:
# allout: (optional) A flag that indicates whether or not to print lines where
#                    a sample is homozygous for the non-archaic allele
#                    (default: don't print non-archaic homozygous sites)
# phased: (optional) The input genotypes are phased, so output one line per
#                    haplotype rather than one line per individual
#                    (default: one line per individual)
# header: (optional) A flag that indicates whether or not to print the header
#                    line (default: don't print header line)
#The "header" line facilitates concatenating the output across chromosomes and
# populations.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   tidprefix="";
   if (length(spop) > 0) {
      tidprefix=spop"_";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", allout);   
   sub(/^[Nn][Oo]?$/, "0", allout);
   if (length(allout) > 0 && allout > 0) {
      print "Printing homozygous non-archaic sites too" > "/dev/stderr";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", phased);   
   sub(/^[Nn][Oo]?$/, "0", phased);
   if (length(phased) > 0 && phased > 0) {
      print "Treating individuals as phased, outputting one line per haplotype variant" > "/dev/stderr";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", header);   
   sub(/^[Nn][Oo]?$/, "0", header);
   if (length(header) > 0 && header > 0) {
      if (length(phased) > 0 && phased > 0) {
         print "CHROM:POS", "Haplotype", "TractID", "TractState";
      } else {
         print "CHROM:POS", "Individual", "TractID", "TractState";
      };
   };
}
#This is just a simple way to enable processing more than 2 distinct input
# files:
FNR==1{
   filenum++;
}
#The first input file is the S' .score file, so we store the column names:
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      sprimecols[$i]=i;
   };
}
#Store the S' tracts in a hash using the archaic alleles as key and the
# tract ID as the value:
#(Note: tract IDs are composed as SPOP_CHROM_SEGMENT or CHROM_SEGMENT if default)
#(Note: ALLELE in the S' .score file is a 0-based index of the allele,
#       analogous to the contents of the GT field in a VCF.)
#We also store the allele index for each position to enable output of
# homozygous non-archaic sites if desired.
filenum==1&&FNR>1{
   loc=$sprimecols["CHROM"] ":" $sprimecols["POS"];
   tid=tidprefix $sprimecols["CHROM"] "_" $sprimecols["SEGMENT"];
   if (length(debug) > 0 && debug >= 1 && tract[loc,tractsites[loc]] != tid) {
      print "Overlapping tracts detected at "loc": "tid" vs. "tract[loc,tractsites[loc]]". Overwriting" > "/dev/stderr";
   };
   tract[loc,$sprimecols["ALLELE"]]=tid;
   tractsites[loc]=$sprimecols["ALLELE"];
}
#The second input file is the output of bcftools query -H -f '%CHROM:%POS[\t%GT]\n'
#Strip the unnecessary components of the bcftools query header, and store
# the column names:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      #The entire header line gets prefixed with "# ":
      sub("# ", "", $i);
      #Each column name is prefixed with it's 1-based index in brackets:
      gsub("[[0-9]+]", "", $i);
      #Each FORMAT field is suffixed with ":[tag]":
      gsub(":GT", "", $i);
      querycols[i]=$i;
   };
}
#Now process the genotypes to identify tract state in each individual at each
# site:
filenum==2&&FNR>1{
   for (i=2; i<=NF; i++) {
      ploidy=split($i, gt, "[/|]", alleledelims);
      #Skip sites not found in the Sprime .score file:
      if (!($1 in tractsites)) {
         continue;
      };
      state="homozygous";
      arcmatch=0;
      tractid=tract[$1,tractsites[$1]];
      #If phased input indicated but site isn't phased, output as
      # state "unphased":
      if (length(phased) > 0 && phased > 0 && alleledelims[1] == "/") {
         for (j=1; j<=ploidy; j++) {
            print $1, querycols[i]"_"j, tractid, "unphased";
         };
         continue;
      };
      for (j=1; j<=ploidy; j++) {
         #If phased input is provided, reset arcmatch for each allele:
         if (length(phased) > 0 && phased > 0) {
            arcmatch=0;
         };
         #The only way for a site to be homozygous is for all alleles to match:
         if (j > 1 && state == "homozygous" && gt[j] != gt[j-1]) {
            state="heterozygous";
         };
         #If any of the alleles match the S' archaic allele, take note:
         if (($1,gt[j]) in tract) {
            arcmatch=1;
            if (length(debug) > 0 && debug >= 1 && tractid != tract[$1,gt[j]]) {
               print "Overlapping tracts detected at "$1": "tractid" vs. "tract[$1,gt[j]] > "/dev/stderr";
            };
            tractid=tract[$1,gt[j]];
         };
         #If phased input was provided, output the tract state for each
         # haplotype:
         if (length(phased) > 0 && phased > 0) {
            #Only output the tract state if it matches the putative archaic
            # allele (per S') unless otherwise requested):
            if (arcmatch) {
               print $1, querycols[i]"_"j, tractid, "archaic";
            } else if (length(allout) > 0 && allout > 0) {
               print $1, querycols[i]"_"j, tractid, "nonarchaic";
            };
         };
      };
      if (length(phased) == 0 || phased < 1) {
         #Only output the tract state if at least one of the alleles matches
         # the putative archaic allele (per S'):
         if (arcmatch) {
            print $1, querycols[i], tractid, state;
         } else if (length(allout) > 0 && allout > 0) {
            print $1, querycols[i], tractid, "homozygous_nonarchaic";
         };
      };
   };
}
