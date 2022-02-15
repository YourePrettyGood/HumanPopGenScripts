#!/bin/awk -f
#This script takes the posterior decoding from admixfrog (the .bin.xz file)
# and performs a very rudimentary run-length encoding that can be used for
# simple visualization of ancestry tracts.
#The output is also a CSV, and has the following columns:
# 1) Chromosome/scaffold ID ("chrom")
# 2) Decoded ancestry state ("state")
# 3) Tract start in SNP space ("id_start")
# 4) Tract end in SNP space ("id_end")
# 5) Tract start in genetic map space ("map_start")
# 6) Tract end in genetic map space ("map_end")
# 7) Tract start in physical/assembly space ("pos_start")
# 8) Tract end in physical/assembly space ("pos_end")
#Starts and ends are specified as 1-based inclusive intervals
# (so like GFF and VCF, not BED or exonerate, etc.)
BEGIN{
   FS=",";
   OFS=FS;
   prevscaf="";
   prevstate="";
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   print "chrom", "state", "id_start", "id_end", "map_start", "map_end", "pos_start", "pos_end";
}
NR>1{
   if ($cols["chrom"] != prevscaf || $cols["viterbi"] != prevstate) {
      if (length(idstart) > 0) {
         print prevscaf, prevstate, idstart, idend, mapstart, mapend, posstart, posend;
      };
      idstart=$cols["id"];
      mapstart=$cols["map"];
      posstart=$cols["pos"];
   };
   prevscaf=$cols["chrom"];
   prevstate=$cols["viterbi"];
   idend=$cols["id"];
   mapend=$cols["map"];
   posend=$cols["pos"];
}
END{
   print prevscaf, prevstate, idstart, idend, mapstart, mapend, posstart, posend;
}
