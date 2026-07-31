#!/bin/awk -f
#This script reformats genetic maps in old SHAPEIT format into the
# format expected by Relate.
#Old SHAPEIT2 format consists of a header line with 4 columns:
# 1) Chromosome
# 2) Physical position (bp)
# 3) Recombination rate (cM/Mbp)
# 4) Genetic position (cM)
#(Ironically, SHAPEIT doesn't actually want the first column...)
#Relate requires genetic map files with a header line with 3 columns:
# 1) Physical position (bp)
# 2) Recombination rate (cM/Mbp)
# 3) Genetic position (cM)
#In addition, Relate genetic map files must be space-separated, not
# tab-separated, whereas the old SHAPEIT files are generally tab-
# separated.
#Optional arguments:
# three_col_input: Input is missing chromosome column
#                  (default=Input has 4 columns, i.e. false or 0)
BEGIN{
   FS="\t";
   OFS=" ";
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", three_col_input);
   sub(/^[Nn][Oo]?$/, "0", three_col_input);
   col_offset=1;
   if (length(three_col_input) > 0 && three_col_input != 0) {
      col_offset=0;
   };
}
NR==1{
   #Print a header line identical to the example at:
   #https://myersgroup.github.io/relate/input_data.html#Data
   print "pos", "COMBINED_rate", "Genetic_Map";
}
NR>1{
   print $(col_offset+1), $(col_offset+2), $(col_offset+3);
}
