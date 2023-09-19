#!/bin/awk -f
#This script takes a series of STDOUT from bcftools stats
# and summarizes overall variant type counts as well as
# singleton counts and counts of substitution classes.
BEGIN{
   #Parse and output tab-separated files:
   FS="\t";
   OFS=FS;
   #Define a map to the substitution classes:
   s["A>C"]="A:T>C:G";
   s["A>G"]="A:T>G:C";
   s["A>T"]="A:T>T:A";
   s["C>A"]="G:C>T:A";
   s["C>G"]="G:C>C:G";
   s["C>T"]="G:C>A:T";
   s["G>A"]="G:C>A:T";
   s["G>C"]="G:C>C:G";
   s["G>T"]="G:C>T:A";
   s["T>A"]="A:T>T:A";
   s["T>C"]="A:T>G:C";
   s["T>G"]="A:T>C:G";
}
#Extract overall variant counts:
/^SN/{
   sub("number of ", "", $3);
   sub(":", "", $3);
   counts[$3]+=$4;
}
#Extract substitution counts:
/^ST/{
   subcounts[s[$3]]+=$4;
   totalsubs+=$4;
}
#Extract singleton counts:
/^SiS/{
   counts["SiSNPs"]+=$4;
   counts["SiIndels"]+=$7;
}
END{
   #Print a header and the overall counts:
   print "Variant Type", "Count", "Percentage";
   print "SNP", counts["SNPs"], counts["SNPs"]*100/counts["records"];
   print "INDELs", counts["indels"], counts["indels"]*100/counts["records"];
   print "multiallelic", counts["multiallelic sites"], counts["multiallelic sites"]*100/counts["records"];
   print "mSNPs", counts["multiallelic SNP sites"], counts["multiallelic SNP sites"]*100/counts["records"];
   print "SiSNPs", counts["SiSNPs"], counts["SiSNPs"]*100/counts["SNPs"];
   print "SiIndels", counts["SiIndels"], counts["SiIndels"]*100/counts["indels"];
   #Print a second header for the substitution class lines:
   print "Substitution Type", "Count", "Percentage";
   #Print the substitution class counts and percentages in order:
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (t in subcounts) {
      print t, subcounts[t], subcounts[t]*100/totalsubs;
   };
}
