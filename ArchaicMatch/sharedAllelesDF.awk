#!/bin/awk -f
#This script outputs a line for each non-archaic, non-ancestral sample
# that match the target population(s) but does not match any excluded
# populations (or are infrequent in the excluded population(s), if
# excludeAF is specified). The output can be further processed in R
# or pre-filtered with subsetVCFstats.pl to identify alleles that
# fall in a BED filter.
# The idea is somewhat generic, but is written from the perspective
# of identifying alleles in modern humans introgressed from archaic
# hominins while mitigating the effects of incomplete lineage sorting
# -- this is the reasoning behind the excludeAF feature.
#As a simplifying assumption to include more sites, the missinghomref
# option is provided to coerce any missing genotypes in the query
# and excluded populations to homozygous reference. This can have a
# fairly large effect on datasets that are merged but not jointly
# genotyped, leading to large segments of mutually exclusive
# genotypes. An expedient but not ideal assumption to make is that
# because no variants were called in this region for the dataset
# with many missing genotypes, the region is invariant, rather than
# having no data to tell either way.
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
   filenum=0;
   #Output a header line:
   print "#CHROM", "POS", "ID", "REF", "ALT", "ALLELE", "ORIGIN", "SAMPLEID", "REGION";
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
   #According to VCF spec, we can assume that GT is the first FORMAT tag,
   # so no need to verify that.
   #Initialize the AN of the excluded class:
   excludeAN=0;
   #By pre-sorting so that targets come first, then excludes, then the rest,
   # we guarantee that all target and exclude alleles are identified
   # before checking the rest of the populations for allele sharing.
   for (i=1; i<=n_indivs; i++) {
      sampleregion=region[ids[indivorder[i]]];
      split($indivorder[i], sample, ":");
      split(sample[1], gt, "[/|]");
      #Skip the individual if it's a target and the genotype is missing,
      # or if it's any other population and we aren't coercing missing
      # genotypes to homozygous reference.
      #If we are coercing missing genotypes to homozygous reference,
      # only do so if the population is not in the targets.
      if (length(missinghomref) == 0 && gt[1] == ".") {
         continue;
      } else if (length(missinghomref) > 0 && gt[1] == "." && sampleregion in targets) {
         continue;
      } else if (length(missinghomref) > 0 && gt[1] == "." && !(sampleregion in targets)) {
         for (j in gt) {
            gt[j]="0";
         };
      };;
      for (allele in gt) {
         samplealleles[gt[allele]]++;
      };
      split($5, alleles, ",");
      alleles[0]=$4;
      for (a in samplealleles) {
         if (sampleregion in targets) {
            targetalleles[sampleregion,a]=1;
            if (length(debug) > 0) {
               print "Target allele", $1, $2, $4, $5, sampleregion, ids[indivorder[i]], a > "/dev/stderr";
            };
         } else if (sampleregion in excludes) {
            if (!(a in excludealleles) && length(debug) > 0) {
               print "Excluded allele", $1, $2, $4, $5, sampleregion, ids[indivorder[i]], a > "/dev/stderr";
            };
            excludealleles[a]+=samplealleles[a];
            excludeAN+=samplealleles[a];
         } else if (!(a in excludealleles) || (excludealleles[a]/excludeAN <= excludeAF)) {
            #We precondition on the allele not matching any excluded alleles.
            #A short-circuit in case no target alleles were identified:
            if (is_empty(targetalleles) == 1) {
               continue;
            };
            #Given a query sample, identify how many target populations
            # the query matches:
            targetmatch="";
            for (j=1; j<=n_targets; j++) {
               if ((targetarr[j],a) in targetalleles) {
                   if (length(targetmatch) == 0) {
                      targetmatch=targetarr[j];
                   } else {
                      targetmatch="Ambiguous";
                   };
               };
            };
            if (length(debug) > 0 && debug > 1) {
               print "Query match", $1, $2, $4, $5, ids[indivorder[i]], a, targetmatch > "/dev/stderr";
            };
            #Now classify the allele based on the number of matches:
            if (length(targetmatch) == 0) {
               #Modern allele, so skip:
               continue;
            } else {
               class=targetmatch;
            };
            #Output a line with site, allele, and sample information:
            print $1, $2, $3, $4, $5, alleles[a], class, ids[indivorder[i]], sampleregion;
#            persamplecount[ids[indivorder[i]],class]+=samplealleles[a];
            if (length(debug) > 0) {
               print "Query match count", $1, $2, $4, $5, ids[indivorder[i]], class, a, persamplecount[ids[indivorder[i]],class] > "/dev/stderr";
            };
         };
      };
      delete samplealleles;
   };
   delete targetalleles;
   if (n_excludes > 0) {
      delete excludealleles;
   };
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
