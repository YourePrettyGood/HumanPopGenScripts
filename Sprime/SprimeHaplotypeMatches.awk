#!/bin/awk -f
#This script takes the S' .score file and the output of bcftools query -H
# -f '%CHROM:%POS[\t%GT]\n' and outputs the phased (or pseudo-phased)
# haplotypic state of each individual polarized by the S' archaic allele.
#The idea is to be able to use the output to plot presence or absence
# of the S' archaic haplotype within a region.
#Note that this script is a bit different from others in that you pass in
# the output of bcftools query twice -- the first time is for the
# archaic samples so that we can output some extra columns indicating
# the state of the archaics relative to the S' allele, and the second time
# is to polarize the rest of the samples. This is simply a consequence of
# not being able to guarantee sample order in the output of bcftools query.
#Options:
# arconly:(optional) A flag that indicates whether to only print sites found
#                    in the S' score file (default: print all sites in VCF)
# header: (optional) A flag that indicates whether or not to print the header
#                    line (default: don't print header line)
#The "header" line facilitates concatenating the output across chromosomes and
# populations.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", arconly);   
   sub(/^[Nn][Oo]?$/, "0", arconly);
   if (length(arconly) > 0 && arconly > 0) {
      print "Only outputting sites found by S'" > "/dev/stderr";
   };
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", header);   
   sub(/^[Nn][Oo]?$/, "0", header);
   if (length(header) > 0 && header > 0) {
      print "CHROM:POS", "Individual", "Haplotype", "TractID", "ArchaicState", "DenisovanAAC", "NeandertalAAC", "AlleleState";
   };
   #For now, we'll hard-code the map from archaic sample ID to group:
   arcpop["Denisova"]="Denisovan";
   arcpop["AltaiNeandertal"]="Neandertal";
   arcpop["Vindija33.19"]="Neandertal";
   arcpop["Chagyrskaya-Phalanx"]="Neandertal";
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
#(Note: tract IDs are composed as CHROM_SEGMENT)
#(Note: ALLELE in the S' .score file is a 0-based index of the allele,
#       analogous to the contents of the GT field in a VCF.)
filenum==1&&FNR>1{
   tract[$sprimecols["CHROM"]":"$sprimecols["POS"],$sprimecols["ALLELE"]]=$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"];
   if (length(arconly) > 0) {
      sprimesites[$sprimecols["CHROM"]":"$sprimecols["POS"]]=$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"];
   };
}
#The second input file is the output of bcftools query -H -f '%CHROM:%POS[\t%GT]\n'
# only including the archaics (i.e. the Denisovan and 3 Neandertals).
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
#Determine the S' archaic allele count in each of the archaic groups:
filenum==2&&FNR>1{
   for (i=2; i<=NF; i++) {
      ploidy=split($i, gt, "[/|]");
      for (j=1; j<=ploidy; j++) {
         if (($1,gt[j]) in tract) {
            arcstate=1;
         } else {
            arcstate=0;
         };
         arcAC[arcpop[querycols[i]]]+=arcstate;
         if (gt[j] != ".") {
            arcAN[arcpop[querycols[i]]]++;
         };
      };
   };
   for (pop in arcAN) {
      if (arcAN[pop] > 0) {
         arcstates[$1,pop]=arcAC[pop];
      } else {
         arcstates[$1,pop]="NA";
      };
   };
   delete arcAC;
   delete arcAN;
}
#The third file is the output of bcftools query -H -f '%CHROM:%POS[\t%GT]\n',
# but this time all the samples are included.
#Strip the unnecessary components of the bcftools query header, and store
# the column names:
filenum==3&&FNR==1{
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
#Now process the phased (or pseudo-phased) genotypes and output each
# haplotypic state for each individual:
filenum==3&&FNR>1{
   #Skip sites not found in S' .score file if requested:
   if (!($1 in sprimesites) && length(arconly) > 0) {
      next;
   };
   for (i=2; i<=NF; i++) {
      ploidy=split($i, gt, "[/|]");
      tractid="";
      #Polarize the haplotypes and identify the tract:
      for (j=1; j<=ploidy; j++) {
         if (($1,gt[j]) in tract) {
            tractid=tract[$1,gt[j]];
            hapstates[j]=1;
         } else if (gt[j] == ".") {
            hapstates[j]="NA";
         } else {
            hapstates[j]=0;
         };
      };
      #In order to have all the same tract ID, we need a separate printing loop:
      for (j=1; j<=ploidy; j++) {
         print $1, querycols[i], j, tractid, hapstates[j], arcstates[$1,"Denisovan"], arcstates[$1,"Neandertal"], gt[j];
      };
      #This delete may not be necessary, since we constrain the loop using ploidy:
#      delete hapstates;
   };
}
