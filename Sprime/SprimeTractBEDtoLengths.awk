#!/bin/awk -f
#This script uses the list of non-outgroup sample IDs plus an arbitrary
# number of tract BED files produced by SprimeTractBED.awk in order to
# output the length of the tract in each individual, plus counts of
# polarized genotypes (only for genotypes containing the archaic allele).
#The output can then be combined with the matchrates file to filter
# tracts and sum up total archaic sequence per individual.
#With the addition of the phased argument, we can do the same per haplotype.
#Arguments:
# pop:    (required) Prefix for the tract ID to make it unique across Sprime
#         runs
# addpop: (optional) Flag indicating whether to add pop as a prefix to tract
#         IDs -- default: Add population prefix (for backwards compatibility)
# phased: (optional) Flag indicating whether or not the input was from
#         phased genotypes
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
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", addpop);   
   sub(/^[Nn][Oo]?$/, "0", addpop);
   if (length(addpop) == 0) {
      addpop=1;
   };
   #ARCMOD is the count of sites heterozygous for the archaic allele
   #ARCARC is the count of sites homozygous for the archaic allele
   #ARC is the count of sites matching the archaic allele when phased
   #MOD is the count of sites not matching the archaic allele when phased
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", header);   
   sub(/^[Nn][Oo]?$/, "0", header);
   sub(/^[Ff]([Aa][Ll][Ss][Ee])?$/, "0", phased);   
   sub(/^[Nn][Oo]?$/, "0", phased);
   if (length(header) > 0 && header > 0) {
      if (length(phased) > 0 && phased > 0) {
         print "CHROM", "Population", "TractID", "Haplotype", "TRACTLEN", "MOD", "ARC";
      } else {
         print "CHROM", "Population", "TractID", "Sample", "TRACTLEN", "ARCMOD", "ARCARC";
      };
   };
   #We hard-code ploidy for haplotype IDs as diploid for now:
   ploidy=2;
}
FNR==NR{
   #Record the IDs of individuals or haplotypes for output ordering:
   if (length(phased) > 0 && phased > 0) {
      for (j=1; j<=ploidy; j++) {
         ids[$1"_"j]=++nindiv;
      };
   } else {
      ids[$1]=++nindiv;
   };
}
FNR<NR{
   split($4, tags, ";");
   for (t in tags) {
      split(tags[t], elems, "=");
      if (elems[1] == "TractID") {
         if (addpop > 0) {
            tractid=pop"_"elems[2];
         } else {
            tractid=elems[2];
         };
      } else if (elems[1] == "Individual" || elems[1] == "Haplotype") {
         id=elems[2];
      } else if (elems[1] == "HomSprimeSites") {
         ArcArc=elems[2];
      } else if (elems[1] == "HetSprimeSites") {
         ArcMod=elems[2];
      } else if (elems[1] == "ArchaicSprimeSites") {
         Arc=elems[2];
      } else if (elems[1] == "ModernSprimeSites") {
         Mod=elems[2];
      };
   };
   if (!(tractid in tractchrom)) {
      tractchrom[tractid]=$1;
   };
   tractlen[tractid,id]=$3-$2;
   tracthet[tractid,id]=ArcMod;
   tracthom[tractid,id]=ArcArc;
   if (length(phased) > 0) {
      tractarc[tractid,id]=Arc;
      tractmod[tractid,id]=Mod;
   };
   if (!(tractid in tracts)) {
      tracts[tractid]=++ntract;
   };
}
END{
   PROCINFO["sorted_in"]="@val_num_asc";
   for (tid in tracts) {
      for (id in ids) {
         if ((tid, id) in tractlen) {
            if (length(phased) > 0 && phased > 0) {
               print tractchrom[tid], pop, tid, id, tractlen[tid,id], tractmod[tid,id], tractarc[tid,id];
            } else {
               print tractchrom[tid], pop, tid, id, tractlen[tid,id], tracthet[tid,id], tracthom[tid,id];
            };
         } else {
            print tractchrom[tid], pop, tid, id, 0, 0, 0;
         };
      };
   };
}
