#!/bin/awk -f

BEGIN{
   FS="\t";
   OFS=FS;
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", phased);
   sub(/^[Nn][Oo]?$/, "0", phased);
}
#CHROM:POS\tIndividual\tTractID\tTractState
#The input is the output of SprimePerSampleTracts.awk, which is a TSV of
# site position, individual ID, S' tract ID, and the genotypic state of
# the tract in that individual. By default, the output of
# SprimePerSampleTracts.awk only contains lines where at least one allele
# matches the S' archaic allele, so this script would only output archaic
# tracts.
#If phased was set for SprimePerSampleTracts.awk, then Individual is actually
# Haplotype, but the code here basically stays the same. The only change is
# that we record the number of archaic vs. non-archaic sites instead of
# homozygous archaic vs. heterozygous sites.
#Be sure to set phased for this script too in that case.
#I've slightly modified this script on 2023/02/21 to set the tract bounds
# as the outer-most positions at which an archaic allele is detected (i.e.
# the start and end are only updated if the current position contains an
# archaic allele). That way the output tracts aren't just fixed to the S'
# tract boundaries if "allout=1" is specified for SprimePerSampleTracts.awk.
!/^CHROM/{
   #tractorder keeps track of the input order of tract IDs so that the output
   # is always consistent with the input:
   if (!($3 in tractorder)) {
      tractorder[$3]=length(tractorder)+1;
      if (length(debug) > 0 && debug > 1) {
         print "New Tract", length(tractorder), $3 > "/dev/stderr";
      };
   };
   #indivorder keeps track of the input order of individuals so that the output
   # is always consistent with the input:
   if (!($2 in indivorder)) {
      indivorder[$2]=length(indivorder)+1;
      if (length(debug) > 0 && debug > 1) {
         print "New Individual", length(indivorder), $2 > "/dev/stderr";
      };
   };
   #Here we start processing each record output by SprimePerSampleTracts.awk:
   split($1, poskey, ":");
   chrom=poskey[1];
   pos=poskey[2];
   #Keep track of the length of the nonarchaic prefixes and suffixes of each
   # archaic haplotype so that we can subtract them off the total site counts.
   # For example, consider a haploid tract NNNAAAANAAANNNNN
   #                                       ppp^^^^^^^^sssss
   # p are part of the prefix, s are part of the suffix, ^ indicates the haplotype
   # As a diploid example: NNNAAAAAAHHNHHNNNNN
   #                       ppp^^^^^^^^^^^sssss
   #The prefix doesn't actually get subtracted, since we reset the count of modern
   # sites to 0 upon encountering the first archaic site.
   #If the tract hasn't been seen before, add it to tracts with initial values:
   if (!(($3,"start",$2) in tracts)) {
      #Only start the tract at a site containing an archaic allele:
      if ($4 != "homozygous_nonarchaic" && $4 != "nonarchaic" && $4 != "unphased") {
         if (length(debug) > 0 && debug > 1) {
            print "Initializing tract "$3" for individual "$2" at "$1" in state "$4 > "/dev/stderr";
         };
         tracts[$3,"chrom",$2]=chrom;
         tracts[$3,"start",$2]=pos;
         tracts[$3,"end",$2]=pos;
         tracts[$3,"state",$2]=$4;
         tracts[$3,"homozygous",$2]=0;
         tracts[$3,"heterozygous",$2]=0;
         tracts[$3,"homozygous_nonarchaic",$2]=0;
         tracts[$3,"archaic",$2]=0;
         tracts[$3,"nonarchaic",$2]=0;
      } else {
         tracts[$3,"prefix",$2]++;
      };
   #If the tract has been seen before, extend it and detect if the state
   # changed:
   } else if ($4 != "unphased") {
      #Move the end of the tract if the current site contains an archaic allele:
      if (pos > tracts[$3,"end",$2]) {
         if ($4 != "homozygous_nonarchaic" && $4 != "nonarchaic") {
            if (length(debug) > 0 && debug > 1) {
               print "Extending tract "$3" for individual "$2" to "pos > "/dev/stderr";
            };
            tracts[$3,"end",$2]=pos;
            #Reset the suffix to 0 every time an archaic allele is found:
            tracts[$3,"suffix",$2]=0;
         } else {
            #Extend the candidate suffix if no archaic allele is found:
            tracts[$3,"suffix",$2]++;
         };
      };
      #Note cases where recombination events likely occurred:
      if ($4 != tracts[$3,"state",$2] && tracts[$3,"state",$2] != "mixed") {
         if (length(debug) > 0) {
            print "Tract "$3" in individual "$2" had state "tracts[$3,"state",$2]" but is now "$4" at position "pos", so coercing to mixed" > "/dev/stderr";
         };
         tracts[$3,"state",$2]="mixed";
      };
   };
   #As a diagnostic for the "mixed" state, we keep track of the number of
   # sites in the tract that the individual had in each state:
   tracts[$3,$4,$2]++;
}
END{
   #Since we put indices as values in the tractorder and indivorder hashes,
   # we set traversal order to "@val_num_asc" to recapitulate input orders:
   PROCINFO["sorted_in"]="@val_num_asc";
   for (t in tractorder) {
      for (i in indivorder) {
         if ((t,"start",i) in tracts) {
            #We're outputting something like BED6 format where the Name
            # column (#4) contains GFF3-like tags:
            if (length(phased) > 0 && phased > 0) {
               print tracts[t,"chrom",i], tracts[t,"start",i]-1, tracts[t,"end",i], "TractID="t";Haplotype="i";State="tracts[t,"state",i]";ArchaicSprimeSites="tracts[t,"archaic",i]";ModernSprimeSites="tracts[t,"nonarchaic",i]-tracts[t,"suffix",i], ".", ".";
            } else {
               print tracts[t,"chrom",i], tracts[t,"start",i]-1, tracts[t,"end",i], "TractID="t";Individual="i";State="tracts[t,"state",i]";HomSprimeSites="tracts[t,"homozygous",i]";HetSprimeSites="tracts[t,"heterozygous",i]";HomModernSites="tracts[t,"homozygous_nonarchaic",i]-tracts[t,"suffix",i], ".", ".";
            };
         };
      };
   };
}
