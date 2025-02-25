#!/bin/awk -f
#This script takes in a tab-separated metadata file and generates the JSON
# file expected by hmmix of ingroup and outgroup sample IDs.
#Outgroup individuals are retained wholesale while ingroup individuals
# can be subsampled sequentially (i.e. not randomly) based on the value(s)
# of target_size.
#Arguments:
# samplecol:   Name of the sample ID column in the input file
# popcol:      Name of the group/pop ID column in the input file
# targets:     Comma-separated list of target/ingroup groups/pops
# outgroups:   Comma-separated list of outgroup groups/pops
# target_size: Comma-separated list of subsample sizes from target groups
#              or a single subsample size reused for each target group
#              (optional, will use all indivs if omitted)
BEGIN{
   FS="\t";
   OFS=FS;
   #Parse the list of outgroups:
   n_outgroups=split(outgroups, OG, ",");
   for (i=1; i<=n_outgroups; i++) {
      OGs[OG[i]]=i;
   };
   #Parse the list of target groups (can be length 1):
   n_targets=split(targets, TG, ",");
   for (i=1; i<=n_targets; i++) {
      TGs[TG[i]]=i;
      #Keep track of the number of indivs kept from each target if
      #target_size is set:
      TGcount[TG[i]]=0;
   };
   #Parse target_size if provided:
   if (length(target_size) > 0) {
      n_sizes=split(target_size, sizes, ",");
      if (n_sizes == 1) {
         #If only one size provided, all targets get the same sample size:
         for (i=1; i<=n_targets; i++) {
            target_sizes[TG[i]]=target_size;
         };
      } else if (n_sizes == n_targets) {
         #If enough sizes are provided, each target gets it's respective sample size:
         for (i=1; i<=n_sizes; i++) {
            target_sizes[TG[i]]=sizes[i];
         };
      } else {
         print "target_size list is not length 1 or the same length as targets, quitting." > "/dev/stderr";
         exit 2;
      };
   };
   #Check if samplecol and popcol are set:
   if (length(samplecol) == 0) {
      print "samplecol (name of sample ID column) not provided, quitting." > "/dev/stderr";
      exit 3;
   };
   if (length(popcol) == 0) {
      print "popcol (name of group ID column) not provided, quitting." > "/dev/stderr";
      exit 4;
   };
}
#Check the header line for column names and make sure samplecol
# and popcol exist:
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!(samplecol in cols)) {
      print "samplecol "samplecol" not found in columns of input, quitting." > "/dev/stderr";
      exit 5;
   };
   if (!(popcol in cols)) {
      print "popcol "popcol" not found in columns of input, quitting." > "/dev/stderr";
      exit 6;
   };
}
#Parse the samples and retain only the requested ones:
NR>1{
   if ($cols[popcol] in TGs) {
      if (length(target_size) == 0 || (length(target_size) > 0 && TGcount[$cols[popcol]] < target_sizes[$cols[popcol]])) {
         ingroup=ingroup",\n    \""$cols[samplecol]"\"";
         TGcount[$cols[popcol]]+=1;
      };
   };
   if ($cols[popcol] in OGs) {
      outgroup=outgroup",\n    \""$cols[samplecol]"\"";
   };
}
END{
   #Generate the JSON:
   printf "{\n  \"ingroup\": [\n%s],\n", substr(ingroup, 3);
   printf "  \"outgroup\": [\n%s]\n}\n", substr(outgroup, 3);
}
