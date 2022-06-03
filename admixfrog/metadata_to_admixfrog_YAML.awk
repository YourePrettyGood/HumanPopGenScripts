#!/bin/awk -f
#This script takes a metadata file and an indication of the parental
# population samples and outgroup sample, and generates an appropriate
# YAML configuration file for Ben Peter's admixfrog.
#
#Input arguments:
# idcol:         (required) Name of the sample ID column in the metadata file
# metacol:       (required) Name of the population ID column in the metadata
#                file
# neandertal:    (optional) Comma-separated list of Neandertal sample IDs
#                (default: AltaiNeandertal,Vindija33.19,Chagyrskaya-Phalanx)
# denisovan:     (optional) Comma-separated list of Denisovan sample IDs
#                (default: Denisova)
# outgroup:      (optional) Comma-separated list of outgroup sample IDs
#                (default: pan_troglodytes)
# pseudohaploid: (optional) Comma-separated list of samples to treat as
#                pseudohaploid
#
#Notes:
# 1) If no outgroup is specified, the default outgroup is automatically
#    added as pseudohaploid
# 2) Certain population IDs are automatically recognized and converted
#    to a short code, any other population IDs are skipped. The
#    recognized populations are:
#     Africa -> AFR
#     America -> AMR
#     CentralAsiaSiberia -> CAS
#     Denisovan -> DEN
#     EastAsia -> EAS
#     Neandertal -> NEA
#     Oceania_other -> OCE
#     Oceania_PIB -> PIB
#     Primates -> ANC
#     SouthAsia -> SAS
#     SoutheastAsia -> SEA
#     WestEurasia -> EUR
#Thus, you'll probably want to modify the code to include your
# particular population IDs and short codes.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(pseudohaploid) == 0) {
      n_pseudohaploid=0;
   } else {
      n_pseudohaploid=split(pseudohaploid, pseudohaploidarr, ",");
      for (i=1; i<=n_pseudohaploid; i++) {
         pseudohaploids[pseudohaploidarr[i]]=i;
      };
   };
   if (length(neandertal) == 0) {
      neandertal="AltaiNeandertal,Vindija33.19,Chagyrskaya-Phalanx";
   };
   n_neandertal=split(neandertal, neandertals, ",");
   if (length(denisovan) == 0) {
      denisovan="Denisova";
   };
   if (length(outgroup) == 0) {
      outgroup="pan_troglodytes";
      if (!(outgroup in pseudohaploids)) {
         pseudohaploids[outgroup]=++n_pseudohaploid;
         pseudohaploidarr[n_pseudohaploid]=outgroup;
      };
   };
   if (length(idcol) == 0) {
      print "Missing idcol variable, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(metacol) == 0) {
      print "Missing metacol variable, please specify it." > "/dev/stderr";
      exit 3;
   };
   #For now, we hard-code the map from long region name to 3-letter code:
   regioncode["Africa"]="AFR";
   regioncode["America"]="AMR";
   regioncode["CentralAsiaSiberia"]="CAS";
   regioncode["Denisovan"]="DEN";
   regioncode["EastAsia"]="EAS";
   regioncode["Neandertal"]="NEA";
   regioncode["Oceania_other"]="OCE";
   regioncode["Oceania_PIB"]="PIB";
   regioncode["Primates"]="ANC";
   regioncode["SouthAsia"]="SAS";
   regioncode["SoutheastAsia"]="SEA";
   regioncode["WestEurasia"]="EUR";
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   region=$cols[metacol];
   #For now we skip any individuals in regions not in the code map
   if (region in regioncode) {
      code=regioncode[region];
      samples[code,++regionsize[code]]=$cols[idcol];
      #Filter down the pseudohaploids provided as a variable based on what's
      # in the metadata file:
      if ($cols[idcol] in pseudohaploids) {
         retained_pseudohaploids[++retained]=$cols[idcol];
      };
   } else {
      print "Skipping sample "$cols[idcol]" because its region ("$cols[metacol]") is not in the region-to-code map." > "/dev/stderr";
   };
}
END{
   print "sampleset:";
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (code in regionsize) {
      printf "    %s:\n", code;
      for (i=1; i<=regionsize[code]; i++) {
         printf "        - %s\n", samples[code,i];
      };
   };
   printf "\n";
   printf "pseudo_haploid:\n";
   for (i=1; i<=retained; i++) {
      printf "    - %s\n", retained_pseudohaploids[i];
   };
}
