#!/bin/awk -f
#This script takes a FASTA index of the reference and the posterior decoding
# from admixfrog (the .bin.xz file) and performs a very rudimentary
# run-length encoding that can be used for simple visualization of ancestry
# tracts.
#The output is a series of BED-format records of tracts in physical/reference
# assembly space, with decoded ancestry state and individual ID in a structured
# Name column (column 4 of standard BED).
#Unlike simpleRLE.awk, this script properly indicates tract ends in physical
# space either based on the next window's start (minus 1), or based on the
# inferred window size and the scaffold lengths from the FASTA index.
#If you set output_type to CSV, it has the following columns instead:
# 1) Chromosome/scaffold ID ("chrom")
# 2) Decoded ancestry state ("state")
# 3) Tract start in SNP space ("id_start")
# 4) Tract end in SNP space ("id_end")
# 5) Tract start in genetic map space ("map_start")
# 6) Tract end in genetic map space ("map_end")
# 7) Tract start in physical/assembly space ("pos_start")
# 8) Tract end in physical/assembly space ("pos_end")
# 9) Sample ID ("sample")
#Starts and ends in this case are specified as 1-based inclusive intervals
# (so like GFF and VCF, not BED or exonerate, etc.)
BEGIN{
   FS=",";
   prevscaf="";
   prevstate="";
   prevpos="";
   if (length(ID) == 0) {
      print "Missing required input: ID" > "/dev/stderr";
      exit 2;
   };
   if (length(output_type) == 0) {
      output_type="BED";
   };
   if (output_type == "CSV") {
      OFS=FS;
   } else {
      OFS="\t";
   };
   #Allow input of the fixed bin size in bp if known.
   #If not input, we set it to 0 until it can be inferred from the data:
   if (length(bin_size) == 0) {
      bin_size=0;
   };
   bin_size_warning=0;
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1{
   n=split($0, a, "\t");
   scaflen[a[1]]=a[2];
}
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (output_type == "CSV") {
      print "chrom", "state", "id_start", "id_end", "map_start", "map_end", "pos_start", "pos_end", "sample";
   };
}
filenum==2&&FNR>1{
   #Infer bin size:
   if (length(prevpos) > 0) {
      if (bin_size == 0) {
         bin_size=$cols["pos"]-prevpos;
      } else if ($cols["pos"]-prevpos != bin_size && bin_size_warning == 0) {
         print "Bin size appears to vary, so scaffold end-fill will use the input or first-inferred size: "bin_size > "/dev/stderr";
         bin_size_warning++;
      };
   };
   #If a tract just terminated, either by end of scaffold or end of state run, output it:
   if ($cols["chrom"] != prevscaf || $cols["viterbi"] != prevstate) {
      #If a tract is terminated by end of state run, 
      if ($cols["chrom"] == prevscaf) {
         posend=$cols["pos"]-1;
      } else {
         #If a tract is terminated by end of scaffold, set tract end to either end of bin or end of scaffold,
         # whichever is smaller:
         if (prevpos+bin_size > scaflen[prevscaf]) {
            posend=scaflen[prevscaf];
         } else {
            posend=prevpos+bin_size-1;
         };
      };
      #Make sure we don't output before the very first line of the input:
      if (length(idstart) > 0) {
         if (output_type == "CSV") {
            print prevscaf, prevstate, idstart, idend, mapstart, mapend, posstart, posend, ID;
         } else {
            print prevscaf, posstart-1, posend, "ID="ID";AncState="prevstate;
         };
      };
      idstart=$cols["id"];
      mapstart=$cols["map"];
      posstart=$cols["pos"];
   };
   prevscaf=$cols["chrom"];
   prevstate=$cols["viterbi"];
   idend=$cols["id"];
   mapend=$cols["map"];
   prevpos=$cols["pos"];
}
END{
   #If a tract is terminated by end of scaffold, set tract end to either end of bin or end of scaffold,
   # whichever is smaller:
   if (prevpos+bin_size > scaflen[prevscaf]) {
      posend=scaflen[prevscaf];
   } else {
      posend=prevpos+bin_size-1;
   };
   if (output_type == "CSV") {
      print prevscaf, prevstate, idstart, idend, mapstart, mapend, posstart, posend, ID;
   } else {
      print prevscaf, posstart-1, posend, "ID="ID";AncState="prevstate;
   };
}
