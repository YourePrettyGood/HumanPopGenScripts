#!/bin/awk -f
#This script is for concatenating IBDNe results from multiple populations
# into a single file based on populations extracted from the input
# filenames.
#The parsing of population from filename is custom, so you may need
# to change it to suit your needs.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
}
FNR==1{
   filenum++;
   #Extract the population from filename:
   n_path_parts=split(FILENAME, path_parts, "/");
   n_fn_parts=split(path_parts[n_path_parts], fn_parts, "_");
   pop=fn_parts[2];
   #Output the header if this is the first file:
   if (filenum == 1) {
      print "Population", "Generation", "Ne", "Lower95PctCI", "Upper95PctCI";
   };
   #Store the column names:
   delete cols;
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
FNR>1{
   print pop, $cols["GEN"], $cols["NE"], $cols["LWR-95%CI"], $cols["UPR-95%CI"];
}
