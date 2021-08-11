#!/bin/awk -f
#This script takes the S' .score file and the output of bcftools query -H
# -f '%CHROM:%POS[\t%GT]\n' and outputs sites in each sample where the
# sample has at least one archaic allele (indicating the state of the
# sample as heterozygous or homozygous).
#Options:
# allout: (optional) A flag that indicates whether or not to print lines where
#                    a sample is homozygous for the non-archaic allele
#                    (default: don't print non-archaic homozygous sites)
# header: (optional) A flag that indicates whether or not to print the header
#                    line (default: don't print header line)
#The "header" line facilitates concatenating the output across chromosomes and
# populations.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   if (length(allout) > 0) {
      print "Printing homozygous non-archaic sites too" > "/dev/stderr";
   };
   if (length(header) > 0) {
      print "CHROM:POS", "Individual", "TractID", "TractState";
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
#(Note: tract IDs are composed as CHROM_SEGMENT)
#(Note: ALLELE in the S' .score file is a 0-based index of the allele,
#       analogous to the contents of the GT field in a VCF.)
filenum==1&&FNR>1{
   tract[$sprimecols["CHROM"]":"$sprimecols["POS"],$sprimecols["ALLELE"]]=$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"];
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
      ploidy=split($i, gt, "[/|]");
      state="homozygous";
      arcmatch=0;
      tractid="";
      #Keep track if the first allele matches the S' archaic allele:
      if (($1,gt[1]) in tract) {
         arcmatch=1;
         tractid=tract[$1,gt[1]];
      };
      for (j=2; j<=ploidy; j++) {
         #The only way for a site to be homozygous is for all alleles to match:
         if (state == "homozygous" && gt[j] != gt[1]) {
            state="heterozygous";
         };
         #If any of the other alleles match the S' archaic allele, take note:
         if (($1,gt[j]) in tract) {
            arcmatch=1;
            tractid=tract[$1,gt[j]];
         };
      };
      #Only output the tract state if at least one of the alleles matches
      # the putative archaic allele (per S'):
      if (arcmatch) {
         print $1, querycols[i], tractid, state;
      } else if (length(allout) > 0) {
         print $1, querycols[i], "NA", "homozygous_nonarchaic";
      };
   };
}
