#!/bin/awk -f
#This script takes in the output of one or more bcftools query calls
# that must include a header and at least the following columns:
# INFO/ID (listed by bcftools query as just ID)
# REF
# ALT
# at least one INFO/AC_* (listed by bcftools query as just AC_*)
#
#Something like the output of:
# bcftools query -H -f "%INFO/ID\\t%INFO/AC_EUR\\t%INFO/AC_ASH\\t%INFO/AC_EAS\\t%INFO/AC_AMR\\t%INFO/AC_AFR\\t%INFO/AC_SAS\\t%INFO/AC_OCN\\t%REF\\t%ALT\\n" HPRCy1_HGSVC3_PIB_mc_CHM13_vcfbub100000_PGprep_biallelic_ACANAF.vcf.gz
#The idea is to output a TSV with one line per variant per population
# it is found in. The same variant may be listed on multiple lines if
# it is shared amongst populations. An extra column is provided with
# a basic classification of the variant type, i.e.:
# SNP: The two alleles differ by a single substitution, but not in length
# INS: The ALT allele is longer than the REF allele and anchor bases match
# DEL: The REF allele is longer than the ALT allele and anchor bases match
# cINS: The ALT allele is longer than the REF, but anchor bases don't match
#       (so there may be a SNP at the anchor base or something more complex)
# cDEL: The REF allele is longer than the ALT, but anchor bases don't match
#       (so there may be a SNP at the anchor base or something more complex)
# MNP: The two alleles differ by multiple substitutions, but not in length
# maybeINV: The two alleles match in length and have less substitutions in
#           the opposite orientation than in the same orientation
#These variant type classifications are not exhaustive or terribly precise.
#Further columns including ALT allele length, allele length difference,
# and the number of substitutions (if SNP, MNP, or maybeINV, NA otherwise)
# in each orientation are also provided.
#
#Input arguments:
# SVthresh: Minimum allele length or length difference to qualify as an SV
#           (default: 50, must be > 0 or will automatically be set to 50)
#           (unit: bp)
BEGIN{
   FS="\t";
   OFS=FS;
   #How long must a variant be in order to qualify as a structural variant?:
   if (length(SVthresh) == 0 || SVthresh <= 0) {
      SVthresh=50;
   };
   print "Using a structural variant length threshold of "SVthresh" bp" > "/dev/stderr";
   n_pops=0;
   PROC_INFO["sorted_in"]="@val_num_asc";
   #Set up a basic hash for performing the complement:
   comp["A"]="T";
   comp["C"]="G";
   comp["G"]="C";
   comp["T"]="A";
   comp["N"]="N";
}
FNR==1{
   for (i=1; i<=NF; i++) {
      #The entire header line gets prefixed with "# ":
      sub("#[ ]?", "", $i);
      #Each column name is prefixed with it's 1-based index in brackets:
      gsub("[[0-9]+]", "", $i);
      #No FORMAT fields are being included, so no need to deal with ":[tag]" suffixes
      querycols[$i]=i;
      #Keep track of the allele count columns specifically:
      if ($i ~ "AC_") {
         n_pops++;
         sub("AC_", "", $i);
         pop_cols[$i]=i;
      };
   };
   #Make sure the minimum four required columns are present:
   if (!("ID" in querycols)) {
      print "ID column missing from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (!("REF" in querycols)) {
      print "REF column missing from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 3;
   };
   if (!("ALT" in querycols)) {
      print "ALT column missing from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 4;
   };
   if (n_pops < 1) {
      print "No AC_* columns detected, cannot proceed. Quitting." > "/dev/stderr";
      exit 5;
   };
   #Print the header line:
   if (FNR==NR) {
      print "ID", "Population", "VariantType", "AlleleLenDiff", "AltAlleleLen", "NumForwardSubstitutions", "NumRevCompSubstitutions";
   };
}
FNR>1{
   #Classify the variant:
   n_alts=split($querycols["ALT"], alts, ",");
   if (n_alts > 1) {
      print "Multiallelic variant detected, input should be decomposed into biallelic variants. Quitting." > "/dev/stderr";
      exit 6;
   };
   isSV=0;
   n_diffs="NA";
   n_rcdiffs="NA";
   #ALT allele length and the difference between ALT and REF allele lengths are useful
   # for classifying variants as SVs or not:
   allelelen=length($querycols["ALT"]);
   SVlen=allelelen-length($querycols["REF"]);
   #Comparing anchor bases between REF and ALT alleles distinguishes simple from complex
   # insertions and deletions:
   anchor_base=substr($querycols["REF"], 1, 1);
   alt_anchor_base=substr($querycols["ALT"], 1, 1);
   if (SVlen == 0) {
      #Start with variants that don't change genome size:
      if (allelelen >= SVthresh) {
         isSV=1;
         #In case the alleles contain substantial flanking sequence, let's do a base-by-base
         # comparison of the two alleles to distinguish SNPs from MNPs and possibly detect
         # very simple inversion candidates:
         split(toupper($querycols["REF"]), refbases, "");
         split(toupper($querycols["ALT"]), altbases, "");
         n_diffs=0;
         n_rcdiffs=0;
         for (j=1; j<=allelelen; j++) {
            if (altbases[j] != refbases[j]) {
               n_diffs+=1;
            };
            #This assumes the two inversion alleles have the same length, which is
            # only likely for young inversions with minimal repeats at breakpoints.
            if (comp[altbases[allelelen-j]] != refbases[j]) {
               n_rcdiffs+=1;
            };
         };
         if (n_diffs == 1) {
            vartype="SNP";
         } else if (n_rcdiffs < n_diffs) {
            vartype="maybeINV";
         } else {
            vartype="MNP";
         };
      };
   } else if (SVlen < 0) {
      #Handle deletions:
      if (-1*SVlen >= SVthresh) {
         isSV=1;
         if (anchor_base == alt_anchor_base) {
            vartype="DEL";
         } else {
            vartype="cDEL";
         };
      };
   } else {
      #Handle insertions:
      if (SVlen >= SVthresh) {
         isSV=1;
         if (anchor_base == alt_anchor_base) {
            vartype="INS";
         } else {
            vartype="cINS";
         };
      };
   };
   if (isSV > 0) {
      #Output a line per structural variant per population it is found in:
      for (pop in pop_cols) {
         if ($pop_cols[pop] > 0) {
            print $querycols["ID"], pop, vartype, SVlen, allelelen, n_diffs, n_rcdiffs;
         };
      };
   };
}
