#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(geno) == 0) {
      print "Missing geno variable specifying the output .geno file, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(snp) == 0) {
      print "Missing snp variable specifying the output .snp file, please specify it." > "/dev/stderr";
      exit 3;
   };
   #Recombination rate scaling factor for physical distance column:
   #Specified as the inverse of rec rate, in bp per cM
   #By default, we use 1M=1000000, since a typical approximation for humans
   # is 1 cM/Mbp
   if (length(recratescale) == 0) {
      recratescale=1000000;
   };
}
NR%2==1{
   sub(":", "_", $1);
   $3=$4/recratescale;
   print $0 >> snp;
}
NR%2==0{
   gsub("0[/|]0", "2");
   gsub("0[/|]1", "1");
   gsub("1[/|]0", "1");
   gsub("1[/|]1", "0");
   gsub("./.", "9");
   gsub("\t", "");
   print $0 >> geno;
}
