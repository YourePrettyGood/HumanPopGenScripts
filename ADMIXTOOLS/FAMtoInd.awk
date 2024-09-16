#!/bin/awk -f
#This script converts PLINK .fam with population in the Family ID column
# into EIGENSTRAT .ind.
BEGIN{
   OFS=" ";
   #Hard-coded column indices for necessary .fam columns:
   famidcol=2;
   fampopcol=1;
   famsexcol=5;
}
#EIGENSTRAT .ind format is space-separated with three columns:
# 1) Sample ID
# 2) Sex encoded as M(=male), F(=female), or U(=unknown)
# 3) Population
{
   sex="U";
   if ($famsexcol == "1") {
      sex="M";
   } else if ($famsexcol == "2") {
      sex="F";
   };
   print $famidcol, sex, $fampopcol;
}
