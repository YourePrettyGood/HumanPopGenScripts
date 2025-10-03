#!/bin/awk -f
#This script generates a simple sample metadata file for stdpopsim
# msprime simulations based on the output of (1) tskit nodes and
# (2) tskit populations on the TreeSequence file output by stdpopsim.
#The sample metadata file has a header with two columns:
# 1) Sample ID ("SampleID")
# 2) Population ID ("Population")
#The latter corresponds to the population name attribute in the
# JSON string output by tskit populations.
#I don't do anything fancy for parsing here, so it may not be robust to
# changes in the format of the output of tskit populations.
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of which file we're on:
   filenum=0;
   #Output a header line:
   print "SampleID", "Population";
}
#Keep track of which file we're on:
FNR==1{
   filenum+=1;
}
#Read in the population indices and names from the table:
filenum==1&&FNR>1{
   #Get rid of the JSON wrapping braces:
   gsub(/[{}]/, "", $2);
   #Parse JSON array elements apart by ", ":
   n=split($2, a, ", ");
   for (i=1; i<=n; i++) {
      #Split the tag-value pair delimited by ": ":
      m=split(a[i], b, ": ");
      #Pull the name out without flanking quotes:
      if (b[1] ~ "name") {
         pop[$1]=substr(b[2], 2, length(b[2])-2);
      };
   };
}
#Parse through the node list, selecting the first node per
# individual and pulling the individual ID and population
# index out, ignoring all unsampled nodes (where population
# == -1):
filenum==2&&FNR>1&&$5>=0{
   if (!($5 in seen)) {
      #Convert individual ID to what is output by tskit vcf, and
      # map population index to population name:
      print "tsk_"$5, pop[$4];
      #Make sure we don't double-print due to the second node for
      # each diploid individual:
      seen[$5]+=1;
   };
}
