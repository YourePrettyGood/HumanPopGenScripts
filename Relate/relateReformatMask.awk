#!/bin/awk -f
#This script converts FASTA-based mask files into the mask format
# expected by Relate.
#Relate mask format is FASTA (may be wrapped) where the only
# permissible sequence characters are P (for passing/included
# sites) and N (for failing/excluded sites).
#Optional arguments:
# pass_char:   Character(s) to consider as passing
#              (default: P)
# trim_header: Trim FASTA header lines of anything after the first space
#              (default: Do not trim)
BEGIN{
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", trim_header);
   sub(/^[Nn][Oo]?$/, "0", trim_header);
   if (length(pass_char) == 0) {
      pass_char="P";
   };
   if (length(trim_header) > 0 && trim_header != 0) {
      
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
