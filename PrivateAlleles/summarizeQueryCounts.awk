#!/bin/awk -f
#This script is somewhat hard-coded, but aggregates the number of
# alleles in certain intersections of VCFs by variant type across
# chromosomes.
#The primary intent is to calculate the number of variants in a
# query callset that are present in databases like dbSNP, the
# 1000 Genomes Project phase 3 callset, and gnomAD.
#I mainly ran this with intersections using dbSNP build 154,
# the 1000 Genomes Project phase 3 callset from 2013/05/02,
# gnomAD exomes from r2.1.1, and gnomAD genomes from r2.1.1.
#This can probably be generalized, but I haven't spent time
# thinking about how to do so.
#The input is generally the output of isecAlleleCounts.awk.
#The input is assumed to have the following columns:
# 1) Chromosome
# 2) Variant subset label
# 3) Variant type
# 4+) Intersection (column name is a binary string indicating which files were
#     intersected)
#
#For my purposes, the second column is either "ALL" to indicate all variants
# or a code indicating variants private to a geographic region (e.g. AFR
# for African-specific variants, EUR for European-specific variants,
# OCN for Oceanian-specific variants).
#Required arguments:
#  subset: Which subset of variants to include (exact match to column 2)
BEGIN{
   FS="\t";
   OFS=FS;
   #This helps iterate over columns in a consistent order:
   PROCINFO["sorted_in"]="@ind_num_asc";
   if (length(subset) == 0) {
      print "Missing subset, please set it. Quitting." > "/dev/stderr";
      exit 2;
   };
}
#For now we assume a particular set of 4 databases used in the intersection,
# so the format of the non-metadata column names is a binary string with 5
# bits indicating the query, dbSNP, 1000 Genomes Project phase 3, gnomAD
# exomes, and gnomAD genomes, respectively.
#Thus, the regexes are hard-coded to account for the intersections involving
# the query and each of these databases.
NR==1{
   for (i=4; i<=NF; i++) {
      #Store the column indices corresponding to intersections including
      # the query:
      if ($i ~ /^1/) {
         qcols[i]=1;
         #Store the subset including both query and dbSNP:
         if ($i ~ /^11/) {
            qdcols[i]=1;
         };
         #Store the subset including both query and 1000 Genomes Project phase 3:
         if ($i ~ /^1[01]1/) {
            qtcols[i]=1;
         };
         #Store the subset including both query and gnomAD exomes:
         if ($i ~ /^1[01][01]1/) {
            qgecols[i]=1;
         };
         #Store the subset including both query and gnomAD genomes:
         if ($i ~ /^1[01][01][01]1/) {
            qggcols[i]=1;
         };
      };
   };
}
NR>1&&$2==subset{
   #Add allele counts from all intersections including query:
   for (i in qcols) {
      q[$3]+=$i;
   };
   #Add allele counts from all intersections including query and dbSNP:
   for (i in qdcols) {
      qd[$3]+=$i;
   };
   #Add allele counts from all intersections including query and 1000GP phase 3:
   for (i in qtcols) {
      qt[$3]+=$i;
   };
   #Add allele counts from all intersections including query and gnomAD exomes:
   for (i in qgecols) {
      qge[$3]+=$i;
   };
   #Add allele counts from all intersections including query and gnomAD genomes:
   for (i in qggcols) {
      qgg[$3]+=$i;
   };
}
END{
   #Sort variant types in lexicographical order:
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Output a line for each variant type:
   for (i in q) {
      print i, q[i], qd[i], qt[i], qge[i], qgg[i];
   };
}
