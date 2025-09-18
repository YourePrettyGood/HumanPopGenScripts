#!/bin/awk -f
#This script takes in a VCF annotated with population-specific allele counts
# and a population of interest and outputs the size of all variants with
# a flag indicating the population-specific variants.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(pop) == 0) {
      print "Missing pop argument, please specify." > "/dev/stderr";
      exit 2;
   };
   ACtag="AC";
   ACpoptag="AC_"pop;
   ANpoptag="AN_"pop;
   #Output the header line:
   print "CHROM", "POS", "VEPALLELE", "TYPE", "ALLELELEN", "SVLEN", "AF_"pop, ACpoptag, pop"_specific";
}
/^#/{
   #print;
}
!/^#/{
   poppresent=0;
   n_info=split($8, infotags, ";");
   for (i=1; i<=n_info; i++) {
      split(infotags[i], tagval, "=");
      if (tagval[1] == ACtag) {
         n_alt=split(tagval[2], AC, ",");
      } else if (tagval[1] == ACpoptag) {
         n_altpop=split(tagval[2], ACpop, ",");
         poppresent+=1;
      } else if (tagval[1] == ANpoptag) {
         ANpop=tagval[2];
         poppresent+=1;
      };
   };
   if (poppresent < 2) {
      print "VCF record is missing needed INFO tags "ANpoptag" and "ACpoptag". Is there a typo in pop="pop"?" > "/dev/stderr";
      print $0 > "/dev/stderr";
      exit 3;
   };
   if (ANpop > 0 && n_alt == n_altpop) {
      for (i=1; i<=n_alt; i++) {
         if (AC[i] == ACpop[i]) {
            if (ACpop[i] > 0) {
               popspecificalts[i]+=1;
            };
         };
      };
   };
   n_alts=split($5, alts, ",");
   if (n_alts != n_alt) {
      print "Mismatch between number of ALT alleles in ALT and number of ALT alleles in "ACtag"." > "/dev/stderr";
      exit 4;
   };
   split($4, refbases, "");
   for (i=1; i<=n_alts; i++) {
      SVlen=length(alts[i])-length($4);
      allelelen=length(alts[i]);
      if (SVlen == 0) {
         split(alts[i], altbases, "");
         n_diffs=0;
         for (j=1; j<=allelelen; j++) {
            if (altbases[j] != refbases[j]) {
               n_diffs+=1;
            };
         };
         if (n_diffs == 1) {
            vartype="SNP";
         } else {
            vartype="MNP";
         };
      } else if (SVlen < 0) {
         if (refbases[1] == substr(alts[i], 1, 1)) {
            vartype="DEL";
         } else {
            vartype="cDEL";
         }
      } else {
         if (refbases[1] == substr(alts[i], 1, 1)) {
            vartype="INS";
         } else {
            vartype="cINS";
         }
      };
      vepallele=alts[i];
      if (length($4) > 1 && vartype !~ "(MNP|cDEL|cINS)") {
         vepallele=substr(vepallele, 2);
         if (vepallele == "") {
            vepallele="-";
         };
      };
      if (ANpop > 0) {
         AF=ACpop[i]/ANpop;
      } else {
         AF="NA";
      };
      if (i in popspecificalts) {
         popspecific="TRUE";
      } else {
         popspecific="FALSE";
      };
      print $1, $2, vepallele, vartype, allelelen, SVlen, AF, ACpop[i], popspecific;
   };
   delete popspecificalts;
}
