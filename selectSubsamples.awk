#!/bin/awk -f
#This script allows you to select samples from a metadata file by either
# randomly subsampling n samples or selecting based on source.
#Note: If you use both filters, be aware that "select" will take precedence,
# so any overlap between the filters may lead to subsamples larger than
# requested.
#
#Options:
# idcol: Column name that contains sample IDs (required)
# samplecol: Column name that contains groups for subsampling (optional)
# seed: PRNG seed to set for reservoir sampling (optional, default=42)
# select: Group(s) to select ("select" means to keep all individuals) (optional)
# selectcol: Column name that contains groups for selecting (optional)
# sizes: Size(s) of subsample(s), can be one size for all, or custom (optional)
# subsample: Group(s) to subsample (optional)
# subselectsize: Subsample size from the selected group(s) (optional)
#
#selectcol is required if select is used
#samplecol and sizes are required if subsample is used
#Using a combination of select and subsample allows for subsampling of some
# groups and keeping all samples from others.
#Using subselectsize requires select, selectcol, and samplecol
#subselectsize leads to subsampling of the select group(s)
#
#Reservoir sampling function based on the pseudocode for Algorithm L at:
# https://en.wikipedia.org/wiki/Reservoir_sampling
#function reservoirsample(s, r, n, k,  i, w, subi) {
#   #Short circuit for n <= k:
#   if (n <= k) {
#      return;
#   };
#   #Due to the way awk instantiates arrays, we will do this step outside the
#   # function:
#   ##Set the initial reservoir:
#   #for (i=1; i<=k; i++) {
#   #   r[i]=s[i];
#   #};
#   #Draw the first random variate to be used for calculating waiting time:
#   w=exp(log(rand())/k);
##   print "Waiting time after item "i" is "w > "/dev/stderr";
#   #Stochastically replace elements of the reservoir until you're through the
#   # input:
#   while (i <= n) {
#      #This step takes advantage of the inverse geometric CDF to move forward
#      # that many items and select the next one as a replacement:
#      #We replace floor() with int(), as S is non-negative.
#      #(S is floor(log(rand())/log(1-w)) )
#      i=i+int(log(rand())/log(1-w))+1;
##      print "Moved to item "i > "/dev/stderr";
#      #Don't fall off the end:
#      if (i <= n) {
#         #Replace a randomly chosen element of the reservoir with the
#         # randomly selected element of the source array:
#         subi=int(k*rand())+1;
#         r[subi]=s[i];
##         print "Replaced item "subi" with item "i" into the reservoir" > "/dev/stderr";
#         #Draw the random variate for calculating the next waiting time:
#         w=exp(log(rand())/k);
##         print "Waiting time after item "i" is "w > "/dev/stderr";
#      };
#   };
#   return;
#}
#Reservoir sampling function based on the pseudocode for Algorithm R at:
# https://en.wikipedia.org/wiki/Reservoir_sampling
function reservoirsample(s, r, n, k,  i, v) {
   #Short circuit for n <= k:
   if (n <= k) {
      return;
   };
   #Due to the way awk instantiates arrays, we will do this step outside the
   # function:
   ##Set the initial reservoir:
   #for (i=1; i<=k; i++) {
   #   r[i]=s[i];
   #};
   for (i=k+1; i<=n; i++) {
      #Draw a random integer from 1 to i:
      v=int(i*rand())+1;
      #Place the sample in the reservoir if it fell in the reservoir:
      if (v <= k) {
         r[v]=s[i];
      };
   };
   return;
}
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(idcol) == 0) {
      print "Missing ID column input variable idcol" > "/dev/stderr";
      exit 2;
   };
   if (length(seed) == 0) {
      seed=42;
   };
   srand(seed);
   if (length(select) == 0 && length(subsample) == 0) {
      print "Neither select nor subsample were selected, no samples will be output" > "/dev/stderr";
      exit 3;
   };
   if (length(select) > 0 && length(selectcol) == 0) {
      print "Please provide the name of the column to select on (selectcol), as it is missing." > "/dev/stderr";
      exit 4;
   };
   if (length(subsample) > 0 && length(samplecol) == 0) {
      print "Please provide the column to group by for subsampling (samplecol), as none was provided." > "/dev/stderr";
      exit 5;
   };
   if (length(subsample) > 0 && length(sizes) == 0) {
      print "Please provide the subsample sizes (sizes), as none were provided." > "/dev/stderr";
      exit 6;
   };
   if (length(select) > 0) {
      n_select=split(select, selectarr, ",");
      for (i=1; i<=n_select; i++) {
         keepall[selectarr[i]]=i;
      };
   };
   if (length(subsample) > 0) {
      n_groups=split(subsample, subsamplearr, ",");
      n_sizes=split(sizes, sizearr, ",");
      for (i=1; i<=n_groups; i++) {
         if (n_sizes != n_groups) {
            #Just take the first subsample size as the size for all, since
            # there aren't enough and I don't want to implement recycling:
            samplesize[subsamplearr[i]]=sizearr[1];
         } else {
            #If there are enough sizes for the groups, assign sizes to
            # the corresponding groups:
            samplesize[subsamplearr[i]]=sizearr[i];
         };
      };
   };
   if (length(subselectsize) > 0) {
      #If we want to subsample from the selected groups too, remove
      # them from keepall and put them in keepsome:
      n_subsizes=split(subselectsize, subsizearr, ",");
      for (i=1; i<=n_select; i++) {
         delete keepall[selectarr[i]];
         if (n_subsizes < n_select) {
            keepsome[selectarr[i]]=subsizearr[1];
         } else {
            keepsome[selectarr[i]]=subsizearr[i];
         };
      };
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   if (length(select) > 0 && $cols[selectcol] in keepall) {
      keep[++num_keep]=$cols[idcol];
      next;
   };
   if (length(select) > 0 && $cols[selectcol] in keepsome) {
      subgroup[++subgroupsize[$cols[selectcol]],$cols[selectcol]]=$cols[idcol];
   };
   if (length(subsample) > 0) {
      group[++groupsize[$cols[samplecol]],$cols[samplecol]]=$cols[idcol];
   };
}
END{
   for (i=1; i<=n_select; i++) {
      if (length(subselectsize) > 0) {
         for (j=1; j<=subgroupsize[selectarr[i]]; j++) {
            subgrouparr[j]=subgroup[j,selectarr[i]];
            #Initialize the reservoir here so awk knows to pass by reference:
            if (j <= keepsome[selectarr[i]]) {
               subgroupsample[j]=subgroup[j,selectarr[i]];
            };
         };
         print "Subsampling "keepsome[selectarr[i]]" from "subgroupsize[selectarr[i]]" for "selectarr[i]" (group #"i" of "n_select")" > "/dev/stderr";
         reservoirsample(subgrouparr, subgroupsample, subgroupsize[selectarr[i]], keepsome[selectarr[i]]);
         reservoir_size=keepsome[selectarr[i]] <= subgroupsize[selectarr[i]] ? keepsome[selectarr[i]] : subgroupsize[selectarr[i]];
         for (j=1; j<=reservoir_size; j++) {
            keep[++num_keep]=subgroupsample[j];
         };
         delete subgrouparr;
         delete subgroupsample;
      };
   };
   for (i=1; i<=n_groups; i++) {
      for (j=1; j<=groupsize[subsamplearr[i]]; j++) {
         grouparr[j]=group[j,subsamplearr[i]];
         #Initialize the reservoir here so awk knows to pass by reference:
         if (j <= samplesize[subsamplearr[i]]) {
            groupsample[j]=group[j,subsamplearr[i]];
         };
      };
      print "Subsampling "samplesize[subsamplearr[i]]" from "groupsize[subsamplearr[i]]" for "subsamplearr[i]" (group #"i" of "n_groups")" > "/dev/stderr";
      reservoirsample(grouparr, groupsample, groupsize[subsamplearr[i]], samplesize[subsamplearr[i]]);
      reservoir_size=samplesize[subsamplearr[i]] <= groupsize[subsamplearr[i]] ? samplesize[subsamplearr[i]] : groupsize[subsamplearr[i]];
      for (j=1; j<=reservoir_size; j++) {
         keep[++num_keep]=groupsample[j];
      };
      delete grouparr;
      delete groupsample;
   };
   PROCINFO["sorted_in"]="@ind_num_asc";
   for (i in keep) {
      print keep[i];
   };
}
