#!/bin/awk -f
#This script takes in the *decoded.*.txt outputs from hmmix decode and
# converts any tracts in the Archaic HMM state with sufficiently high
# posterior probability (i.e. >= minposterior) into BED format.
#This is useful for computing overlaps with introgressed tract inferences
# from other methods (e.g. ArchaicSeeker2.0, Sprime projections, S*)
# using something like BEDtools.
#Arguments:
# minposterior: Minimum threshold for posterior probability of an Archaic tract
# hapid:        Haplotype ID of the input file
# nea_refs:     Comma-separated list of Neanderthal reference VCF sample IDs
#               Order in this list determines order of columns in the output
#               (default: Vindija33.19,Chagyrskaya-Phalanx,AltaiNeandertal)
# den_refs:     Comma-separated list of Denisovan reference VCF sample IDs
#               Order in this list determines order of columns in the output
#               (default: Denisova)
#2026/10/16: Output sample ID and extra -admixpop columns if they exist
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(minposterior) == 0) {
      print "minposterior variable is missing, please set it." > "/dev/stderr";
      print "minposterior is the minimum threshold of the posterior probability of HMM state for an Archaic tract." > "/dev/stderr";
      exit 2;
   };
   if (minposterior < 0.0 || minposterior > 1.0) {
      print "minposterior "minposterior" is outside of [0,1], please fix it." > "/dev/stderr";
      exit 3;
   };
   if (length(nea_refs) == 0) {
      nea_refs="Vindija33.19,Chagyrskaya-Phalanx,AltaiNeandertal";
   };
   n_nea=split(nea_refs, nea_arr, ",");
   for (i=1; i<=n_nea; i++) {
      nea[nea_arr[i]]=i;
   };
   if (length(den_refs) == 0) {
      den_refs="Denisova";
   };
   n_den=split(den_refs, den_arr, ",");
   for (i=1; i<=n_den; i++) {
      den[den_arr[i]]=i;
   };
   PROCINFO["sorted_in"]="@val_num_asc";
}
#Map column names to indices:
FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   has_nea=1;
   for (n in nea) {
      if (!(n in cols)) {
         has_nea=0;
      };
   };
   has_den=1;
   for (d in den) {
      if (!(d in cols)) {
         has_den=0;
      };
   };
}
FNR>1{
   if ($cols["state"] == "Archaic" && $cols["mean_prob"] >= minposterior) {
      if ("admixpopvariants" in cols) {
         if (has_nea > 0) {
            nea_matches="";
            for (n in nea) {
               nea_matches=nea_matches""OFS""$cols[n];
            };
            if (has_den > 0) {
               den_matches="";
               for (d in den) {
                  den_matches=den_matches""OFS""$cols[d];
               };
               print $cols["chrom"], $cols["start"], $cols["end"], hapid, $cols["mean_prob"], $cols["admixpopvariants"], substr(nea_matches, 2), substr(den_matches, 2);
            } else {
               print $cols["chrom"], $cols["start"], $cols["end"], hapid, $cols["mean_prob"], $cols["admixpopvariants"], substr(nea_matches, 2);
            };
         } else {
            print $cols["chrom"], $cols["start"], $cols["end"], hapid, $cols["mean_prob"];
         };
      } else {
         print $cols["chrom"], $cols["start"], $cols["end"], hapid, $cols["mean_prob"];
      };
   };
}
