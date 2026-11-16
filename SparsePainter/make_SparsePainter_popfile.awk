#!/bin/awk -f
#This script takes in a sample metadata file in tab-separated format (with header)
# and a file simply listing sample IDs in the order found in a VCF (e.g. the output
# of `bcftools query -l` on the VCF), and outputs a tab-separated file (no header)
# with two columns: sample ID, and population label.
#This is the format expected by SparsePainter for the -popfile argument.
#If the `all_vs_all` argument is set to 1, the population label column will be
# set to the population label followed by an underscore and an integer, producing
# a unique population label for every sample.
#Options:
# samplecol:  Name of the sample ID column in the metadata file
#             (default: SampleID)
# popcol:     Name of the population label column in the metadata file
# all_vs_all: Whether to augment the provided population labels with an integer
#             to produce sample-unique labels required for all-vs-all mode of
#             SparsePainter
#             (default: use provided population labels, i.e. target-vs-ref or ref-vs-ref modes)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(samplecol) == 0) {
      print "No samplecol provided. Defaulting to SampleID." > "/dev/stderr";
      samplecol="SampleID";
   };
   if (length(popcol) == 0) {
      print "Please provide the name of the population column in the metadata file. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(all_vs_all) == 0) {
      all_vs_all=0;
   };
   #Keep track of which file we're on:
   filenum=0;
}
FNR==1{
   filenum+=1;
}
#First file is the metadata file:
filenum==1&&FNR==1{
   #Check the header for necessary columns:
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!("SampleID" in cols && popcol in cols)) {
      print "Metadata file "FILE" is missing one or more of the expected columns: "samplecol" and "popcol". Please check that you input the correct file. Quitting." > "/dev/stderr";
      exit 3;
   };
}
filenum==1&&FNR>1{
   popmap[$cols[samplecol]]=$cols[popcol];
}
#Second file is the list of samples in the order of the target VCF:
filenum==2{
   if (!($1 in popmap)) {
      print "Sample "$1" is missing from the metadata file. Quitting." > "/dev/stderr";
      exit 4;
   };
   if (all_vs_all > 0) {
      popn[popmap[$1]]+=1;
      poplabel=popmap[$1]"_"popn[popmap[$1]];
   } else {
      poplabel=popmap[$1];
   };
   print $1, poplabel;
}
