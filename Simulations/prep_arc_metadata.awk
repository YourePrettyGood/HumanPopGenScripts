#!/bin/awk -f
#This script takes any number of sample metadata files produced by
# ts_to_metadata.awk and does some simple reformatting for
# compatibility with archaicMatchSprime.awk and SprimeArchaicMatchRate.awk.
#The two main reformatting items are:
# 1) Rename the columns to "Sample" and "Region"
# 2) Extract rows with Population matching values in simgroups and
#    replace the population label with the corresponding value in
#    groups
#Arguments:
# simgroups: Comma-separated list of archaic species/groups in the
#            simulation sample metadata file(s)
#            (default: NeaA,DenA)
# groups:    Comma-separated list of archaic species/group names to
#            use for compatibility with downstream match rate scripts
#            Order must match the order of simgroups, as they will
#            be used as 1:1 replacements
#            (default: Neandertal,Denisovan)
BEGIN{
   FS="\t";
   OFS=FS;
   #Get the lists of groups to rename and do some validation:
   if (length(simgroups) == 0) {
      simgroups="NeaA,DenA";
      print "Using default simgroups, which corresponds to PapuansOutOfAfrica_10J19" > "/dev/stderr";
   };
   n_simgroups=split(simgroups, simgrouparr, ",");
   if (length(groups) == 0) {
      groups="Neandertal,Denisovan";
      print "Using default groups: Neandertal,Denisovan" > "/dev/stderr";
   };
   n_groups=split(groups, grouparr, ",");
   if (n_simgroups != n_groups) {
      print "Length of simgroups ("n_simgroups") does not match length of groups ("n_groups"). Cannot proceed." > "/dev/stderr";
      exit 2;
   };
   for (i=1; i<=n_simgroups; i++) {
      group_map[simgrouparr[i]]=grouparr[i];
   };
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum+=1;
   #Also validate that the input files have appropriate headers as expected:
   if ($1 != "SampleID" || $2 != "Population") {
      print "Header for file "filenum" ("FILENAME") header is missing an expected column. Cannot proceed." > "/dev/stderr";
      exit 3;
   };
}
#Read in each file, keeping count of each sample and extracting those
# with a Population found in the group_map:
FNR>1{
   sample_count[$1]+=1;
   if ($2 in group_map) {
      if ($1 in arc_samples && arc_samples[$1] != group_map[$2]) {
         print "Sample "$1" ("arc_samples[$1]") found in another metadata file in a different Population ("group_map[$2]"). Cannot proceed." > "/dev/stderr";
         exit 4;
      };
      arc_samples[$1]=group_map[$2];
   };
}
#Now validate the inputs all matched and output the archaic metadata file:
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (id in sample_count) {
      if (sample_count[id] != filenum) {
         print "Sample "id" not found in all input files. Cannot proceed." > "/dev/stderr";
         exit 5;
      };
   };
   #Output the archaic metadata file:
   print "Sample", "Region";
   for (id in arc_samples) {
      print id, arc_samples[id];
   };
}
