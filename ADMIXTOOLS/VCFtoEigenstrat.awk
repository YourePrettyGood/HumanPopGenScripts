#!/bin/awk -f
#This script takes in a genetic map (or multiple per-chromosome genetic maps)
# as well as the output of bcftools query -H -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
# and produces .snp and .geno files for EIGENSTRAT or similar programs
# (e.g. ASCEND, ADMIXTOOLS, etc.).
#If you don't have a genetic map available, you can alternatively specify
# the value of "recrate" in cM/Mbp as a constant recombination rate scaling
# factor. The genetic position of each site will thus be computed as:
# g = x * 1000000 / r
# where r is the value of recrate and x is the value of the POS column.
#This script understands two main formats of genetic maps:
# 1) shapeit (has a header and three columns: pos, chr, and cM)
# 2) plink (aka beagle, simply the PLINK .map format)
#          (no header, four columns: CHROM, ID, CM, and POS)
#These correspond to two of the main formats for human genetic maps,
# i.e. those provided for SHAPEIT4 and Beagle 5.4
#The headers get ignored, so make sure the columns are in the appropriate order.
#Mandatory arguments:
# outprefix:    Prefix to use for the output .snp and .geno files
# numgenmaps:   The number of genetic map files passed in
#               (required if genetic map files are passed in)
# genmapformat: A string indicating the format of the genetic map files
#               (one of "shapeit", "plink", or "beagle")
#               (required if genetic map files are passed in)
# recrate:      Fixed recombination rate to use in place of a genetic map
#               (specified in cM/Mbp)
#               (required if *no* genetic map files are passed in)
#Optional arguments:
# setmissingvarids: Flag to indicate whether to set ID if the VCF has it as "."
#                   (default: No, ID format: [CHROM]_[POS]_[REF]_[ALT])
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(outprefix) == 0) {
      print "Missing prefix for output .snp and .geno files, please specify it." > "/dev/stderr";
      exit 2;
   };
   snpfile=outprefix".snp";
   genofile=outprefix".geno";
   #Recombination rate is specified either with genetic map files
   # or a constant recombination rate:
   if (length(numgenmaps) > 0 && numgenmaps > 0) {
      print "Using "numgenmaps" genetic map files passed before the bcftools query output." > "/dev/stderr";
      #Check which format of genetic map we have:
      if (length(genmapformat) == 0 || tolower(genmapformat) !~ /plink|beagle|shapeit/) {
         print "No or invalid genetic map format specified. Please specify plink, beagle, or shapeit." > "/dev/stderr";
         exit 3;
      };
   } else if (length(recrate) > 0 && recrate > 0) {
      print "Scaling physical positions into genetic map positions using fixed recombination rate "recrate" cM/Mbp" > "/dev/stderr";
      recratescale=1000000/recrate;
      numgenmaps=0;
   } else {
      print "No valid parameters specified to determine genetic map positions. Please either specify numgenmaps and pass that many genetic maps before the bcftools query output, or else specify recrate in cM/Mbp and do not pass any genetic maps." > "/dev/stderr";
      exit 4;
   };
   #Do we want to set variant IDs when they aren't set?:
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", setmissingvarids);
   sub(/^[Nn][Oo]?$/, "0", setmissingvarids);
   if (length(setmissingvarids) > 0 && setmissingvarids > 0) {
      print "Setting missing variant IDs using pattern [CHROM]_[POS]_[REF]_[ALT]" > "/dev/stderr";
   };
   #Initialize the map_index to 0, so the maps use 1-based indexing:
   map_index=0;
   #Also initialize prev_index, which keeps track of the last map marker we
   # used:
   prev_index=0;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
  filenum++;
}
filenum<=numgenmaps{
   #Because awk arrays are associative arrays, we have to do a bit of trickery
   # to be able to efficiently search for the closest genetic map markers:
   #We actually use four arrays:
   # 1) Index to physical position (pos_map)
   # 2) Chromosome to start index/offset (chrom_map)
   # 3) Index to genetic position (cM_map)
   # 4) Chromosome to end index/offset (end_map)
   #We assume that the input genetic map is never out of order and contiguous
   # within a chromosome, which seems a safe assumption.
   #Increment the map_index, since we're on a new marker:
   map_index+=1;
   if (tolower(genmapformat) == "shapeit") {
      #The genetic map files provided with SHAPEIT4 have a header, so skip:
      if ($1 == "pos") {
         next;
      };
      #The format is %POS\t%CHROM\t%CM:
      #Only set the chrom_map offset if it doesn't already exist:
      if (!($2 in chrom_map)) {
         chrom_map[$2]=map_index;
      };
      #end_map moves along until the last marker for the chromosome:
      end_map[$2]=map_index;
      #Set the values of pos_map and cM_map:
      pos_map[map_index]=$1;
      cM_map[map_index]=$3;
   } else if (tolower(genmapformat) ~ /plink|beagle/) {
      #The genetic map files provided with BEAGLE5.4 have no header, and
      # the format is: %CHROM\t%ID\t%CM\t%POS
      #This is PLINK's .map format, so not exactly specific to genetic maps.
      #Only set the chrom_map offset if it doesn't already exist:
      if (!($1 in chrom_map)) {
         chrom_map[$1]=map_index;
      };
      #end_map moves along until the last marker for the chromosome:
      end_map[$1]=map_index;
      #Set the values of pos_map and cM_map:
      pos_map[map_index]=$4;
      cM_map[map_index]=$3;
   };
   #For interpolation, we perform linear interpolation like PLINK1.9 does in
   # the apply_cm_map() function. This means taking the flanking genetic map
   # markers, taking the discrete derivative of genetic position with respect
   # to physical position, scaling this discrete derivative by the physical
   # distance from the closest preceding marker, and adding this to the genetic
   # position of the closest preceding marker.
}
filenum>numgenmaps&&FNR==1{
   n_gts=0;
   for (i=1; i<=NF; i++) {
      sub("#[ ]?", "", $i);
      gsub("[[0-9]+]", "", $i);
      if ($i ~ ":GT") {
         gsub(":GT", "", $i);
         n_gts+=1;
         gtcols[i]=n_gts;
      } else {
         cols[$i]=i;
      };
   };
   if (!(("ID" in cols) && ("CHROM" in cols) && ("POS" in cols) && ("REF" in cols) && ("ALT" in cols))) {
      print "bcftools query output is missing one or more of ID, CHROM, POS, REF, or ALT. Cannot proceed." > "/dev/stderr";
      exit 5;
   };
   if (n_gts == 0) {
      print "No genotype columns found in bcftools query output. Cannot proceed." > "/dev/stderr";
      exit 6;
   };
}
#EIGENSTRAT .geno format uses the count of REF/ANC alleles, and 9 for missing.
filenum>numgenmaps&&FNR>1{
   if (length(numgenmaps) > 0 && numgenmaps > 0) {
      #First we compute the physical position of the current site:
      #Find a marker in the genetic map that is either the same as the current
      # site or the closest marker *after* the current site:
      #Emit a warning and skip the site if it's not on a chromosome in the
      # genetic map:
      if (!($cols["CHROM"] in chrom_map)) {
         print "Site "$cols["CHROM"]":"$cols["POS"]" unable to be processed, as chromosome "$cols["CHROM"]" is not in the provided genetic map. Skipping this site." > "/dev/stderr";
         next;
      };
      #Skip to the correct chromosome:
      if (prev_index < chrom_map[$cols["CHROM"]]) {
         map_index=chrom_map[$cols["CHROM"]];
      };
      #Iterate forward until map marker position is either geq VCF site position, or
      # we've reached the end of the chromosome:
      while (pos_map[map_index] < $cols["POS"] && map_index < end_map[$cols["CHROM"]]) {
         map_index+=1;
      };
      #Yes, three of the following cases have identical actions, but they represent
      # distinct cases, so the code logic is easier to follow this way.
      if (pos_map[map_index] == $cols["POS"]) {
         #Site exactly overlaps a map marker, so no interpolation necessary:
         gen_pos=cM_map[map_index];
      } else if (map_index == chrom_map[$cols["CHROM"]]) {
         #Site is before the first map marker on the chromosome, so no interpolation:
         #(Basically, set gen_pos to 0.0)
         gen_pos=cM_map[map_index];
      } else if (map_index == end_map[$cols["CHROM"]] && $cols["POS"] > pos_map[map_index]) {
         #Site is after the last map marker on the chromosome, so no interpolation:
         #(Set gen_pos to the genetic position of the last map marker)
         gen_pos=cM_map[map_index];
      } else {
         #Any other site, which is thus found between a pair of map markers,
         # so we need to interpolate:
         #Calculate the slope between the map markers:
         dg_dx=(cM_map[map_index] - cM_map[map_index-1])/(pos_map[map_index] - pos_map[map_index-1]);
         #Determine the physical distance from the left map marker to the query site:
         delta_x=($cols["POS"] - pos_map[map_index-1]);
         #Find the genetic position of the left map marker:
         g_0=cM_map[map_index-1];
         #Do the linear interpolation:
         #g_i = g_0 + \delta x \cdot \frac{dg}{dx}
         #where
         #\delta x = x_i - x_0
         #\frac{dg}{dx} = (g_1 - g_0)/(x_1 - x_0)
         #for left map marker 0 and right map marker 1, and query site i 
         gen_pos=g_0 + delta_x*dg_dx;
      };
      if (length(debug) > 0 && debug > 1) {
         print $cols["CHROM"], $cols["POS"], pos_map[map_index-1], pos_map[map_index], map_index, chrom_map[$cols["CHROM"]], end_map[$cols["CHROM"]], prev_index, gen_pos > "/dev/stderr";
      };
      #Set prev_index for the next iteration:
      prev_index=map_index;
   } else if (recratescale > 0) {
      gen_pos=$cols["POS"]*recratescale;
      if (length(debug) > 0 && debug > 1) {
         print $cols["CHROM"], $cols["POS"], gen_pos, recratescale > "/dev/stderr";
      };
   };
   #Add an ID if the VCF didn't have one:
   if (length(setmissingvarids) > 0 && setmissingvarids > 0 && $cols["ID"] == ".") {
      $cols["ID"]=$cols["CHROM"]"_"$cols["POS"]"_"$cols["REF"]"_"$cols["ALT"];
   };
   #Now print the .snp line:
   print $cols["ID"], $cols["CHROM"], gen_pos, $cols["POS"], $cols["REF"], $cols["ALT"] >> snpfile;
   #Compose the .geno line:
   genoline="";
   for (i=1; i<=NF; i++) {
      if (i in gtcols) {
         #Re-initialize gt to 0 for each individual:
         gt=0;
         #Split up the VCF genotype into an array of alleles of with ploidy elements:
         ploidy=split($i, alleles, "[/|]");
         #Add 1 to gt for each REF allele in the genotype,
         # though set gt to 9 if any allele is missing.
         #This means gt may be >9 if the genotype is only partially missing like ./0
         #This only encounters a problem if ploidy is >=9, so most folks should be fine.
         for (j=1; j<=ploidy; j++) {
            if (alleles[j] == ".") {
               gt=9;
            } else if (alleles[j] == "0") {
               gt+=1;
            };
         };
         #Convert gt>9 back to gt==9 to be compatible with EIGENSTRAT .geno format:
         if (gt >= 9) {
            gt=9;
         };
         genoline=genoline""gt;
      };
   };
   #Now print the .geno line:
   print genoline > genofile;
}
