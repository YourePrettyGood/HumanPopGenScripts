#!/bin/awk -f
#This script counts alleles for each non-archaic, non-ancestral sample
# that match the target population(s) but do not match any excluded
# populations (or are infrequent in the excluded population(s), if
# excludeAF is specified). The idea is somewhat generic, but is
# written from the perspective of identifying alleles in modern humans
# introgressed from archaic hominins while mitigating the effects of
# incomplete lineage sorting -- this is the reasoning behind the
# excludeAF feature.
#As a simplifying assumption to include more sites, the missinghomref
# option is provided to coerce any missing genotypes in the query
# and excluded populations to homozygous reference. This can have a
# fairly large effect on datasets that are merged but not jointly
# genotyped, leading to large segments of mutually exclusive
# genotypes. An expedient but not ideal assumption to make is that
# because no variants were called in this region for the dataset
# with many missing genotypes, the region is invariant, rather than
# having no data to tell either way.
#A newer feature is to also filter on whether or not the site has
# an ID (which is a typical way to assess novel vs. known sites).
# Specifying the novelonly flag will trigger this feature.
#
#Options/arguments:
# idcol:        Name of the column in the metadata file with sample IDs
#               (required)
# metacol:      Name of the column in the metadata file with population IDs
#               (required)
# target:       Comma-separated list of population IDs which are the targets
#               of allele sharing (i.e. the output is counts of alleles
#               shared between query sample individuals and the target
#               population)
# exclude:      Comma-separated list of population IDs used to identify
#               alleles that are candidates for ILS, and thus excluding
#               those alleles
# excludeAF:    Maximum allele frequency threshold in the excluded population(s)
#               allowed before omitting an allele
#               (default: omit any allele (0.0))
# missinghomref:Coerce missing genotypes to homozygous reference
#               (only occurs in samples from the query and excluded pops)
#               (default: don't coerce (unset))
# novelonly:    Only assess sites where the ID column of the VCF is "."
#               (default: assess all sites (unset))
#A helper function to perform the popcount instruction, i.e. count the
# number of bits set in an integer:
#Based on: https://rosettacode.org/wiki/Population_count#AWK
function pop_count(xx,  xq, y) {
   while (xx > 0) {
      xq=int(xx/2);
      if (xx-xq*2 == 1) {
         y++;
      };
      xx=xq;
   };
   return y;
}
#Another helper function for finding which bit is selected, given
# an integer with pop_count() == 1:
#Note: 1-based bit indexing so we can easily plug it into an awk array
function which_bit(x,  b, y) {
   y=0;
   b=2^y;
   while (x >= b) {
      if (x == b) {
         return y+1;
      };
      y++;
      b=2^y;
   };
   print "which_bit() was called with an int that wasn't a power of 2" > "/dev/stderr";
}
#A helper function to check if an array is empty:
#Based on: https://stackoverflow.com/questions/20075845/how-to-check-if-awk-array-is-empty
function is_empty(a,  i) {
   for (i in a) {
      return 0;
   };
   return 1;
}
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(idcol) == 0) {
      print "Missing idcol variable, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(metacol) == 0) {
      print "Missing metacol variable, please specify it." > "/dev/stderr";
      exit 3;
   };
   if (length(target) == 0) {
      print "Missing target population for allele sharing, please specify it." > "/dev/stderr";
      exit 4;
   } else {
      print "Targeting population(s) "target" for allele sharing." > "/dev/stderr";
   };
   n_targets=split(target, targetarr, ",");
   n_classes=0;
   for (i in targetarr) {
      targets[targetarr[i]]=i;
      classes[++n_classes]=targetarr[i];
   };
   classes[++n_classes]="Ambiguous";
#   classes[++n_classes]="Modern";
   if (length(exclude) > 0) {
      print "Excluding alleles shared with populations: "exclude > "/dev/stderr";
      n_excludes=split(exclude, excludearr, ",");
      for (i in excludearr) {
         excludes[excludearr[i]]=i;
      };
   } else {
      print "No populations set for exclusion." > "/dev/stderr";
      n_excludes=0;
   };
   if (length(excludeAF) == 0) {
      excludeAF=0.0;
   };
   if (n_excludes > 0) {
      print "Alleles shared with excluded populations will only be retained if their frequency is "excludeAF" or less." > "/dev/stderr";
   };
   if (length(missinghomref) > 0) {
      print "Converting missing genotypes in excluded and query populations to homozygous reference." > "/dev/stderr";
   };
   if (length(novelonly) > 0) {
      print "Only counting allele at novel sites (ID column is .)" > "/dev/stderr";
   };
   filenum=0;
}
FNR==1{
   filenum++;
}
#First file is the metadata file
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
#Keep track of the metacol value for each sample:
filenum==1&&FNR>1{
   region[$cols[idcol]]=$cols[metacol];
}
#Second file is the VCF
#Organize sample iteration order prior to the sample records so that
# the process is simplified:
filenum==2&&/^#CHROM/{
   #We establish two maps: ids and regionmap
   #These maps, plus the regionsize array, allow us to iterate over regions
   # slowly and individuals within a region quickly.
   for (i=10; i<=NF; i++) {
      ids[i]=$i;
#      regionsize[region[$i]]++;
      regionmap[region[$i],++regionsize[region[$i]]]=i;
   };
   #Now determine the order of the regions, starting with target(s) and then
   # exclude(s), then the rest:
   n_indivs=0;
   for (p=1; p<=n_targets; p++) {
      for (i=1; i<=regionsize[targetarr[p]]; i++) {
         indivorder[++n_indivs]=regionmap[targetarr[p],i];
      };
   };
   for (p=1; p<=n_excludes; p++) {
      for (i=1; i<=regionsize[excludearr[p]]; i++) {
         indivorder[++n_indivs]=regionmap[excludearr[p],i];
      };
   };
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (p in regionsize) {
      if (!(p in targets) && !(p in excludes)) {
         for (i=1; i<=regionsize[p]; i++) {
            indivorder[++n_indivs]=regionmap[p,i];
         };
      };
   };
   if (length(debug) > 0 && debug > 2) {
      for (i=1; i<=n_indivs; i++) {
         print indivorder[i], $indivorder[i], region[$indivorder[i]] > "/dev/stderr";
      };
   };
}
#Now, on a per-record basis, identify which samples and regions share alleles
# with the target populations but not the excluded populations
filenum==2&&!/^#/{
   #Skip the site if novelonly is set and ID isn't "." (i.e. missing):
   if (length(novelonly) > 0 && $3 != ".") {
      next;
   };
   #Print a definition line of the site:
   printf ">%s\t%s\t%s\t%s\n", $1, $2, $4, $5;
   sampleregion="";
   n_alts=split($5, alleles, ",");
   alleles[0]=$4;
   #According to VCF spec, we can assume that GT is the first FORMAT tag,
   # so no need to verify that.
   #Initialize the AN of the excluded class:
   excludeAN=0;
   #By pre-sorting so that targets come first, then excludes, then the rest,
   # we guarantee that all target and exclude alleles are identified
   # before checking the rest of the populations for allele sharing.
   for (i=1; i<=n_indivs; i++) {
      #Print a prefix with the region once per region:
      if (region[ids[indivorder[i]]] != sampleregion) {
         if (sampleregion != "") {
            printf "\n";
         };
         printf "%s\t", region[ids[indivorder[i]]];
      };
      sampleregion=region[ids[indivorder[i]]];
      split($indivorder[i], sample, ":");
      ploidy=split(sample[1], gt, "[/|]");
      #Skip the individual if it's a target and the genotype is missing,
      # or if it's any other population and we aren't coercing missing
      # genotypes to homozygous reference.
      #If we are coercing missing genotypes to homozygous reference,
      # only do so if the population is not in the targets.
      if (length(missinghomref) == 0 && gt[1] == ".") {
         printf "%iX", ploidy; #Missing data (X)
         continue;
      } else if (length(missinghomref) > 0 && gt[1] == "." && sampleregion in targets) {
         printf "%iX", ploidy; #Missing data (X)
         continue;
      } else if (length(missinghomref) > 0 && gt[1] == "." && !(sampleregion in targets)) {
         for (j in gt) {
            gt[j]="0";
         };
      };;
      for (allele in gt) {
         samplealleles[gt[allele]]++;
      };
      for (a in samplealleles) {
         if (sampleregion in targets) {
            targetalleles[sampleregion,a]=1;
            if (length(debug) > 0) {
               print "Target allele", $1, $2, $4, $5, sampleregion, ids[indivorder[i]], a > "/dev/stderr";
            };
            printf "%i%s", samplealleles[a], alleles[a];
         } else if (sampleregion in excludes) {
            if (!(a in excludealleles) && length(debug) > 0) {
               print "Excluded allele", $1, $2, $4, $5, sampleregion, ids[indivorder[i]], a > "/dev/stderr";
            };
            excludealleles[a]+=samplealleles[a];
            excludeAN+=samplealleles[a];
            printf "%iE", samplealleles[a]; #Excluded allele
#            #We're excluding these alleles from consideration as shared,
#            # so the fastest solution is to remove them from the target
#            # allele set:
#            for (t in targets) {
#               delete targetalleles[t,a];
#            };
         } else if (!(a in excludealleles) || (excludealleles[a]/excludeAN <= excludeAF)) {
            #We precondition on the allele not matching any excluded alleles.
            #A short-circuit in case no target alleles were identified:
            if (is_empty(targetalleles) == 1) {
               printf "%iY", samplealleles[a]; #EmptY targetalleles
               continue;
            };
            #Given a query sample, identify how many target populations
            # the query matches:
 #           targetmatch=0;
            targetmatch="";
            for (j=1; j<=n_targets; j++) {
               if ((targetarr[j],a) in targetalleles) {
#                  targetmatch+=1*2^(j-1);
                  if (length(targetmatch) == 0) {
                     targetmatch=substr(targetarr[j], 1, 1);
                  } else {
                     targetmatch="A"; #Ambiguous
                  };
               };
            };
            if (length(debug) > 0 && debug > 1) {
               print "Query match", $1, $2, $4, $5, ids[indivorder[i]], a, targetmatch > "/dev/stderr";
            };
            #Now classify the allele based on the number of matches:
#            if (targetmatch == 0) {
#               #Modern allele:
##               class="Modern";
#               continue;
#            } else if (pop_count(targetmatch) == 1) {
#               #Only one target match:
#               class=targetarr[which_bit(targetmatch)];
#            } else {
#               #Ambiguous, matches multiple targets:
#               class="Ambiguous";
#            };
#            persamplecount[ids[indivorder[i]],class]+=samplealleles[a];
            if (length(targetmatch) == 0) {
               printf "%iM", samplealleles[a]; #Modern allele
            } else {
               printf "%i%s", samplealleles[a], targetmatch;
            };
            if (length(debug) > 0) {
               print "Query match count", $1, $2, $4, $5, ids[indivorder[i]], class, a, persamplecount[ids[indivorder[i]],class] > "/dev/stderr";
            };
         } else {
            printf "%iE", samplealleles[a]; #Excluded allele
         };
      };
      delete samplealleles;
   };
   printf "\n";
   delete targetalleles;
   if (n_excludes > 0) {
      delete excludealleles;
   };
#   delete popalleles;
}
#END{
#   for (i=1; i<=n_indivs; i++) {
#      for (c=1; c<=n_classes; c++) {
#         class=classes[c];
#         if (!((ids[indivorder[i]],class) in persamplecount)) {
#            persamplecount[ids[indivorder[i]],class]=0;
#         };
#         print ids[indivorder[i]], region[ids[indivorder[i]]], class, persamplecount[ids[indivorder[i]],class];
#      };
#   };
#}
