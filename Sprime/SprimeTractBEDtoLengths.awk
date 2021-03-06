#!/bin/awk -f
#This script uses the list of non-outgroup sample IDs plus an arbitrary
# number of tract BED files produced by SprimeTractBED.awk in order to
# output the length of the tract in each individual, plus counts of
# polarized genotypes (only for genotypes containing the archaic allele).
#The output can then be combined with the matchrates file to filter
# tracts and sum up total archaic sequence per individual.
#Arguments:
# pop:    (required) Prefix for the tract ID to make it unique across Sprime
#         runs
# header: (optional) Flag indicating whether or not to output the header
#         line
BEGIN{
   FS="\t";
   OFS=FS;
   nindiv=0;
   ntract=0;
   if (length(pop) == 0) {
      print "Missing pop variable, please set it." > "/dev/stderr";
      print "pop is used as a prefix for the tract ID to make it unique" > "/dev/stderr";
      print "across Sprime runs" > "/dev/stderr";
      exit 2;
   };
   #ARCMOD is the count of sites heterozygous for the archaic allele
   #ARCARC is the count of sites homozygous for the archaic allele
   if (length(header) > 0) {
      print "CHROM", "Population", "TractID", "Sample", "TRACTLEN", "ARCMOD", "ARCARC";
   };
}
FNR==NR{
   ids[$1]=++nindiv;
}
FNR<NR{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], elems, "=");
      if (elems[1] == "TractID") {
         tractid=pop"_"elems[2];
      } else if (elems[1] == "Individual") {
         id=elems[2];
      } else if (elems[1] == "HomSprimeSites") {
         ArcArc=elems[2];
      } else if (elems[1] == "HetSprimeSites") {
         ArcMod=elems[2];
      };
   };
   if (!(tractid in tractchrom)) {
      tractchrom[tractid]=$1;
   };
   tractlen[tractid,id]=$3-$2;
   tracthet[tractid,id]=ArcMod;
   tracthom[tractid,id]=ArcArc;
   if (!(tractid in tracts)) {
      tracts[tractid]=++ntract;
   };
}
END{
   PROCINFO["sorted_in"]="@val_num_asc";
   for (tid in tracts) {
      for (id in ids) {
         if ((tid, id) in tractlen) {
            print tractchrom[tid], pop, tid, id, tractlen[tid,id], tracthet[tid,id], tracthom[tid,id];
         } else {
            print tractchrom[tid], pop, tid, id, 0, 0, 0;
         };
      };
   };
}
