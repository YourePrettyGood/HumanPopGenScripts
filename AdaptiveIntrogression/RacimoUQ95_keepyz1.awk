#!/bin/awk -f
#This script takes a VCF of the archaics processed with bcftools +fill-tags
# -- -t AC,AF,AN and identifies sites with the appropriate (y,z) pattern
# for the Racimo et al. 2017 U_{A,B,C,D}(w,x,y,z) and Q95_{A,B,C,D}(x,y,z)
# statistics. That is, given C=NEA and D=DEN, (y,z)==(1.0,1.0) is AMB,
# (y,z)==(0.0,1.0) is DEN, and (y,z)==(1.0,0.0) is NEA, where y is the
# frequency of the derived allele in C and z is the DAF in D.
#For now, we set this up as creating separate files from separate runs
# of this script, one for each classified archaic origin: AMB, DEN, or NEA.
#As such, the mandatory argument "origin" controls which (y,z) threshold is
# applied to filter the input VCF.
#The output is a filtered VCF intended for use in filtering a VCF of modern
# human genotypes in order to calculate w, x, and the U and Q95 statistics.
#The filtered VCF should be compressed with bgzip, indexed with tabix, and
# then can be used with the -T argument to bcftools view to filter the modern
# human VCF for only these sites (although make sure to also include -v snps
# to avoid incidental indels). This approach significantly shortens processing
# time and disk space usage compared to merging archaic and modern VCFs and
# then calculating frequencies, though it requires a bit more care when
# performing the join between archaic and modern frequency outputs.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(origin) == 0) {
      print "Missing origin, must be set to AMB, DEN, or NEA. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (origin !~ /^(AMB|DEN|NEA)$/) {
      print "Invalid origin specified ("origin"), must be set to AMB, DEN, or NEA. Quitting." > "/dev/stderr";
      exit 3;
   };
}
/^#/{
   print;
}
!/^#/{
   #Construct the set of alleles:
   if ($5 != ".") {
      nalt=split($5, alleles, ",");
   };
   alleles[0]=$4;
   nalleles=nalt+1;
   #Initialize the two frequency arrays:
   DENAF["A"]=0.0;
   DENAF["C"]=0.0;
   DENAF["G"]=0.0;
   DENAF["T"]=0.0;
   NEAAF["A"]=0.0;
   NEAAF["C"]=0.0;
   NEAAF["G"]=0.0;
   NEAAF["T"]=0.0;
   #Parse the INFO string:
   ninfo=split($8, info, ";");
   #Make sure we set a sentinel in case AA not found or non-ACGT:
   AA=".";
   #Keep track of whether AF tags are available, as fully homozygous ref
   # sites apparently don't get AF set, as there are no ALT alleles:
   DENAFset=0;
   NEAAFset=0;
   #Find and parse the INFO/AA, INFO/AF_DEN, and INFO/AF_NEA tags:
   for (i in info) {
      split(info[i], tagparts, "=");
      if (tagparts[1] == "AA" && tagparts[2] ~ /[ACGTacgt]/) {
         AA=tagparts[2];
      } else if (tagparts[1] == "AF_DEN") {
         #Store the DEN ALT allele frequencies if found:
         sumdenaf=0.0;
         ndenaf=split(tagparts[2], denafs, ",");
         for (i=1; i<=ndenaf; i++) {
            sumdenaf+=denafs[i];
            DENAF[toupper(alleles[i])]=denafs[i];
         };
         #Calculate REF allele frequency for DEN:
         DENAF[toupper(alleles[0])]=1.0-sumdenaf;
         #Mark that DENAF has been set:
         DENAFset=1;
      } else if (tagparts[1] == "AF_NEA") {
         #Store the NEA ALT allele frequencies if found:
         sumneaaf=0.0;
         nneaaf=split(tagparts[2], neaafs, ",");
         for (i=1; i<=nneaaf; i++) {
            sumneaaf+=neaafs[i];
            NEAAF[toupper(alleles[i])]=neaafs[i];
         };
         #Calculate REF allele frequency for NEA:
         NEAAF[toupper(alleles[0])]=1.0-sumneaaf;
         #Mark that NEAAF has been set:
         NEAAFset=1;
      };
   };
   #Fill in the REF allele frequency if no AF tag was provided:
   if (!DENAFset) {
      DENAF[toupper($4)]=1.0;
   };
   if (!NEAAFset) {
      NEAAF[toupper($4)]=1.0;
   };
   #If the ancestral allele call is usable, output the line:
   if (AA ~ /[ACGTacgt]/) {
      #In particular, only output the line if the ancestral allele frequency
      # in DEN and/or NEA is 0.0:
      #This assumes that the site is biallelic, as the criterion we're actually
      # going for is derived allele frequency of 1.0 in DEN or NEA.
      #We split this out by desired origin:
      # AMB is DEN_AAF == 0.0 && NEA_AAF == 0.0
      # DEN is DEN_AAF == 0.0 && NEA_AAF == 1.0
      # NEA is DEN_AAF == 1.0 && NEA_AAF == 0.0
      DEN_AAF=DENAF[toupper(AA)];
      NEA_AAF=NEAAF[toupper(AA)];
      if (origin == "AMB" && DEN_AAF == 0.0 && NEA_AAF == 0.0) {
         print $0;
      } else if (origin == "DEN" && DEN_AAF == 0.0 && NEA_AAF == 1.0) {
         print $0;
      } else if (origin == "NEA" && DEN_AAF == 1.0 && NEA_AAF == 0.0) {
         print $0;
      };
   };
}
