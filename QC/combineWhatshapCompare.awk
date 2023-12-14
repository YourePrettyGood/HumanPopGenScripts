#!/bin/awk -f
#This script summarizes the phasing error rates calculated by
# `whatshap compare` across chromosomes.
#Some aspects of this script are very customized, like
# extracting the assembly from the filename.
BEGIN{
   FS="\t";
   OFS=FS;
   #Print header line:
   print "Sample asm truth,query", "Total SER", "SER", "FER", "SFER", "numpairs", "Blockwise HER", "covsites";
}
FNR==1{
   delete cols;
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (FILENAME ~ /GRCh38/) {
      asm="GRCh38";
   } else {
      asm="hs37d5";
   };
}
FNR>1{
   #Index the array by sample, assembly, and query,truth:
   i=$cols["#sample"]" "asm" "$cols["dataset_name0"]","$cols["dataset_name1"];
   #
   t[i]+=$cols["all_assessed_pairs"];
   split($cols["all_switchflips"], sf, "/");
   s[i]+=sf[1];
   f[i]+=sf[2];
   c[i]+=$cols["covered_variants"];
   b[i]+=$cols["blockwise_hamming"];
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (i in t) {
      if (t[i] > 0) {
         #Total Switch Error Rate (TSER), i.e. naïvely handling flips
         tser=(s[i]+2*f[i])*100/t[i];
         #Switch Error Rate (SER), i.e. ignoring flips
         ser=s[i]*100/t[i];
         #Flip Error Rate (FER), i.e. ignoring non-flip switches
         fer=f[i]*100/t[i];
         #Switch-Flip Error Rate (SFER), i.e. correctly handling flips
         sfer=(s[i]+f[i])*100/t[i];
         #Blockwise Hamming Error Rate (BHER)
         bwher=b[i]*100/c[i];
         #Print the line:
         print i, tser, ser, fer, sfer, t[i], bwher, c[i];
      };
   };
}
