#!/bin/awk -f
#This script takes a metadata file and generates the .poplabels file
# used by Relate.
#Your metadata file is expected to have a single header line followed
# by sample lines, and should be whitespace-delimited.
#This means that .fam is not directly supported, since it is missing
# a header line.
#I may add support in the future based on specifying column indices
# instead of names.
#.poplabels format has a header line and is space-delimited:
# 1) sample - Sample ID (must match input)
# 2) population - Population ID
# 3) group - Superpopulation/Region ID
# 4) sex - 1 for male, 2 for female, NA for unknown, use 1 if haploid
#Note that Relate uses NA for diploid and 1 if haploid.
#Mandatory arguments:
# idcol:     Column name in the input for the sample ID
# popcol:    Column name in the input for the population ID
# groupcol:  Column name in the input for the group ID
# sexcol:    Column name in the input for the sample sex
#Optional arguments:
# ploidycol: Column name in the input for the sample ploidy
BEGIN{
   OFS=" ";
   if (length(idcol) == 0) {
      print "Missing idcol, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(popcol) == 0) {
      print "Missing popcol, please specify it." > "/dev/stderr";
      exit 3;
   };
   if (length(groupcol) == 0) {
      print "Missing groupcol, please specify it." > "/dev/stderr";
      exit 4;
   };
   if (length(sexcol) == 0) {
      print "Missing sexcol, please specify it." > "/dev/stderr";
      exit 5;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!(idcol in cols)) {
      print "idcol ("idcol") missing from input metadata file." > "/dev/stderr";
      exit 2;
   };
   if (!(popcol in cols)) {
      print "popcol ("popcol") missing from input metadata file." > "/dev/stderr";
      exit 3;
   };
   if (!(groupcol in cols)) {
      print "groupcol ("groupcol") missing from input metadata file." > "/dev/stderr";
      exit 4;
   };
   if (!(sexcol in cols)) {
      print "sexcol ("sexcol") missing from input metadata file." > "/dev/stderr";
      exit 5;
   };
   if (length(ploidycol) > 0 && !(ploidycol in cols)) {
      print "ploidycol ("ploidycol") specified but missing from input metadata file." > "/dev/stderr";
      exit 6;
   };
   print "sample", "population", "group", "sex";
}
NR>1{
   #Convert sex column value to 1 if male, 2 if female, or NA if unknown:
   #Some typical alternate value schemes:
   # PLINK .fam/.psam: 1 = male, 2 = female, 0 = unknown
   # 1kGP .panel: male = male, female = female
   # HGDP metadata: M = male, F = female
   # PIBv1 metadata: M = male, F = female
   # Human genetic sex: XY = male, XX = female
   #Note that this section will fail if your metadata file doesn't follow
   # one of these conventions.
   #We aren't handling sex chromosome aneuploidies here either.
   if ($cols[sexcol] ~ /^(male|M|XY)$/) {
      recoded_sex=1;
   } else if ($cols[sexcol] ~ /^(female|F|XX)$/) {
      recoded_sex=2;
   } else {
      recoded_sex="NA";
   };
   if (length(ploidycol) > 0 && $cols[ploidycol] < 2) {
      recoded_sex=1;
   } else {
      recoded_sex="NA";
   };
   print $cols[idcol], $cols[popcol], $cols[groupcol], recoded_sex;
}
