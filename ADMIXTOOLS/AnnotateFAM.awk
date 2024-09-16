#!/bin/awk -f
#This script annotates a PLINK .fam file with sex and population information
# from a tab-separated metadata file.
#The first input file is the metadata file, the second is the existing .fam
# file.
#Population is placed in the "Family ID" column of the .fam, as is customary
# for some datasets.
#We recognize a limited number of sex encodings in the metadata file.
#The main ones are:
# - M or F
# - male or female
# - Male or Female
# - XY or XX
#Obviously, these are incomplete encodings, but they cover
# several commonly-used cases like 1000 Genomes Project and
# HGDP.
#If neither of the options match, sex is set to 0 (Unknown).
#Mandatory arguments:
# idcol:  Name of the column containing sample ID in the metadata file
# popcol: Name of the column containing population in the metadata file
# sexcol: Name of the column containing sex in the metadata file
BEGIN{
   FS="\t";
   OFS=" ";
   if (length(idcol) == 0) {
      print "Missing idcol variable, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(popcol) == 0) {
      print "Missing popcol variable, please specify it." > "/dev/stderr";
      exit 3;
   };
   if (length(sexcol) == 0) {
      print "Missing sexcol variable, please specify it." > "/dev/stderr";
      exit 4;
   };
   #Hard-coded values for .fam columns of interest:
   famidcol=2;
   fampopcol=1;
   famsexcol=5;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
filenum==1&&FNR==1{
   #Map the metadata file column names to column indices:
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Check that the necessary columns exist in the metadata file:
   if (!(idcol in cols)) {
      print "ID column "idcol" missing from metadata file." > "/dev/stderr";
      exit 2;
   };
   if (!(popcol in cols)) {
      print "Population column "popcol" missing from metadata file." > "/dev/stderr";
      exit 3;
   };
   if (!(sexcol in cols)) {
      print "Sex column "sexcol" missing from metadata file." > "/dev/stderr";
      exit 4;
   };
}
#Store the population and sex of each individual in the metadata file:
filenum==1&&FNR>1{
   pop[$cols[idcol]]=$cols[popcol];
   sex[$cols[idcol]]=$cols[sexcol];
}
#Annotate individuals in the .fam with population and sex:
filenum==2{
   #Handle space-separated .fam if necessary:
   if (NF == 1) {
      n=split($1, a, " ");
      for (i=1; i<=n; i++) {
         $i=a[i];
      };
   };
   #Replace the population and sex if found in metadata:
   if ($famidcol in pop) {
      new_sex="0";
      if (sex[$famidcol] ~ /^([Mm]([Aa][Ll][Ee])?|XY)$/) {
         new_sex="1";
      } else if (sex[$famidcol] ~ /^([Ff]([Ee][Mm][Aa][Ll][Ee])?|XX)$/) {
         new_sex="2";
      };
      $famsexcol=new_sex;
      $fampopcol=pop[$famidcol];
      npopspaces=gsub(" ", "_", $fampopcol);
      if (npopspaces > 0 && !($fampopcol in replacedpops)) {
         print "Replaced spaces in "pop[$famidcol]" with underscores" > "/dev/stderr";
         replacedpops[$fampopcol]+=1;
      };
   };
   print $0;
}
