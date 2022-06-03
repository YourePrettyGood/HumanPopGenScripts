#!/bin/awk -f
#This script merges the CSVs produced by `smcpp plot --csv` across runs
# or iterations, allowing for re-plotting and co-plotting of
# inferred Ne(t) trajectories in R.
#The code for parsing the phasing state and inputs is custom to the
# filenames I was using, but the general idea is that your filenames
# should be underscore-delimited, with the second token being the
# inputs, and the second-to-last token being the phasing state.
#By "inputs", I mean whether a single distinguished lineage was
# used, or the type of composite likelihood used
# (e.g. composite10 would be composite likelihood from 10
#  distinguished lineages).
#By "phasing state", I mean whether the input data were unphased
# or phased.
#SMC++ calls the final model "iteration 0", so this script converts
# iteration 0 to "final".
BEGIN{
   FS=",";
   OFS=FS;
   filenum=0;
   #Print a header:
   print "Population", "Time", "Ne", "Iteration", "Phasing", "Inputs";
}
FNR==1{
   filenum++;
   n_path_parts=split(FILENAME, patharr, "/");
   n_fn_parts=split(patharr[n_path_parts], fnarr, "_");
   inputs=fnarr[2];
   phasing=fnarr[n_fn_parts-2];
}
FNR>1{
   gsub(/\r/, "");
   iteration=$5 == 0 ? "final" : $5;
   print $1, $2, $3, iteration, phasing, inputs;
}
