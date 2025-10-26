#!/bin/awk -f
#This script takes a sample metadata file produced by ts_to_metadata.awk
# and performs two tasks:
# 1) Relabels the Population column for samples from selected populations
#    aside from a subset of selected size
# 2) Adds a Region column corresponding to the original population names or
#    substitute region names specified by the user for the groups subject
#    to subsetting
#As an example, let's say your input metadata file contains 100 samples
# from a population labeled CEU and you want to create a subset of 25
# and add a Region column for them called EUR. Simply specify:
# -v "subsets=CEU" -v "subsetsizes=25" -v "regions=EUR"
# and the result will be that the first 25 CEU samples will retain
# "CEU" in the Population column, the following 75 CEU samples will
# be relabeled as "notCEU" in the Population column, and all 100
# formerly CEU samples will have "EUR" in the Region column.
# Samples from any other population will retain their Population
# label, and their Region column will be the same as their Population
# column.
#I'm mainly using this to address the discrepancy between the
# PapuansOutOfAfrica_10J19 model simulating four extant populations
# (YRI, CEU, CHB, Papuan) while our empirical datasets included
# larger superpopulations/regions (e.g. AFR, EUR, EAS, OCN) as well
# as populations (e.g. Yoruba, French, Han, Goroka), so the simulations
# effectively need to simulate the regions, but we need to be able
# to easily pull out subsets corresponding to what we actually used
# for analyses (e.g. Yoruba as Sprime outgroup, Goroka as Sprime and
# selection scan target, OCN used for LD calculations).
#This is a bit of an oversimplification, as there's definitely
# structure in the empirical data that isn't modeled by the simulations.
#Required arguments:
# subsets:     Comma-separated list of Population names to split based
#              on subset size and rename
# subsetsizes: Comma-separated list of subset sizes in the same order as
#              the populations in "subsets"
# regions:     Comma-separated list of names to use in the new Region
#              column for the split populations in the same order as the
#              populations in "subsets"
BEGIN{
   FS="\t";
   OFS=FS;
   #Validate the input arguments:
   if (length(subsets) == 0) {
      print "List of populations to subset (subset) not specified, cannot proceed." > "/dev/stderr";
      exit 2;
   };
   n_subset_pops=split(subsets, subsetarr, ",");
   for (i=1; i<=n_subset_pops; i++) {
      subsetpops[subsetarr[i]]=i;
   };
   if (length(subsetsizes) == 0) {
      print "List of subset sizes (subsetsizes) not specified, cannot proceed." > "/dev/stderr";
      exit 3;
   };
   n_subset_sizes=split(subsetsizes, subsetsizearr, ",");
   #Double-check that the user didn't mismatch the lists:
   if (n_subset_sizes != n_subset_pops) {
      print "subsets and subsetsizes do not contain the same number of items, cannot proceed." > "/dev/stderr";
      exit 4;
   };
   for (i=1; i<=n_subset_sizes; i++) {
      subsetcount[subsetarr[i]]=subsetsizearr[i];
   };
   if (length(regions) == 0) {
      print "List of region names for subsets and their overall sets (regions) not specified, cannot proceed." > "/dev/stderr";
      exit 5;
   };
   n_subset_regions=split(regions, regionarr, ",");
   #We're more flexible about the regions list, as too many isn't
   # really a problem:
   if (n_subset_regions > n_subset_pops) {
      print "regions and subsets do not contain the same number of items, ignoring extras." > "/dev/stderr";
   } else if (n_subset_regions < n_subset_pops) {
      print "regions does not contain enough items to match subsets, cannot proceed." > "/dev/stderr";
      exit 6;
   };
   #Extra regions in the list are ignored, only the first n_subset_pops
   # are actually used:
   for (i=1; i<=n_subset_pops; i++) {
      regionmap[subsetarr[i]]=regionarr[i];
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Make sure the required input columns exist:
   if (!("SampleID" in cols && "Population" in cols)) {
      print "SampleID or Population column missing from input metadata, cannot proceed." > "/dev/stderr";
      exit 7;
   };
   #Print the new header:
   print "SampleID", "Population", "Region";
}
NR>1{
   if ($cols["Population"] in subsetpops) {
      #Retain the original population label if in the initial subset:
      if (subsetcount[$cols["Population"]] > 0) {
         print $cols["SampleID"], $cols["Population"], regionmap[$cols["Population"]];
         subsetcount[$cols["Population"]] -= 1;
      } else {
      #Otherwise, prefix the population label with "not":
         print $cols["SampleID"], "not"$cols["Population"], regionmap[$cols["Population"]];
      };
   } else {
      #Use the population label as the region label for any populations
      # not subset:
      print $cols["SampleID"], $cols["Population"], $cols["Population"];
   };
}
