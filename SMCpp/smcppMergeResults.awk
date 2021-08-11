#!/bin/awk -f

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
