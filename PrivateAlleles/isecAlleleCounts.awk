#!/bin/awk -f
#This script takes the output of bcftools isec -c none -n+1 and counts
# variants in each intersection in total and by variant type.
#Note that the bcftools isec command does not include -p, as we want
# a 5-column list of variants as input.
#Required arguments:
#  chrom:  Which chromosome the input represents
#  subset: The variant subset the input represents
#          (e.g. ALL for all variants, AFR for African-specific variants)
#  query:  Alias for subset
BEGIN{
   OFS="\t";
   #Check for required arguments:
   if (length(chrom) == 0) {
      print "Missing chrom, please set it. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(query) == 0 && length(subset) == 0) {
      print "Missing query, please set it. Quitting." > "/dev/stderr";
      exit 3;
   };
   if (length(subset) == 0) {
      subset=query;
   };
}
{
   #Get the lengths of the REF and ALT alleles to determine variant type:
   rl=length($3);
   al=length($4);
   #Increment the total variant count for the given intersection:
   c[$5]++;
   #Also increment the by-type variant counts for the given intersection:
   if ($4 ~ /</) {
      #ALT will have a < if it is a symbolic allele
      bytype[$5,"SYMBOLIC"]++;
   } else if (rl == al && rl == 1) {
      #Both REF and ALT are of length 1 if SNP
      bytype[$5,"SNP"]++;
   } else if (rl == al && rl > 1) {
      #Both REF and ALT are the same length and length > 1 if MNP
      bytype[$5,"MNP"]++;
   } else if (rl < al && rl == 1) {
      #INS has length(REF) of 1 and ALT is longer than REF
      bytype[$5,"INS"]++;
   } else if (rl < al && rl > 1) {
      bytype[$5,"COMPLEX"]++;
   } else if (rl > al && al == 1) {
      #DEL has length(ALT) of 1 and REF is longer than ALT
      bytype[$5,"DEL"]++;
   } else {
      bytype[$5,"COMPLEX"]++;
   };
}
END{
   #Iterate over intersections in a consistent order:
   PROCINFO["sorted_in"]="@ind_num_asc";
   #For now, we hard-code the list of variant types:
   split("SNP,INS,DEL,MNP,COMPLEX,SYMBOLIC", varianttype, ",");
   #Print a header line:
   printf "CHROM\tSubset\tVariantType";
   for (i in c) {
      printf "\t%s", i;
   };
   printf "\n";
   #Print the total variant intersection counts:
   printf "%s\t%s\tALL", chrom, subset;
   for (i in c) {
      printf "\t%i", c[i];
   };
   printf "\n";
   for (t in varianttype) {
      #Print the type-specific intersection counts:
      printf "%s\t%s\t%s", chrom, subset, varianttype[t];
      for (i in c) {
         printf "\t%i", bytype[i,varianttype[t]];
      };
      printf "\n";
   };
}
