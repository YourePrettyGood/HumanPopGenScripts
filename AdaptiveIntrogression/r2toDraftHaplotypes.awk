#!/bin/awk -f
#This script takes the output of vcftools --hap-r2 (with --min-r2 set)
# on a single scaffold/chromosome and outputs a list of the input SNPs
# annotated with the cluster they belong to and the number of edges
# missing from the graph of that cluster (where edges represent LD
# levels passing the --min-r2 filter).
#The idea is to replicate the first clustering step of Gittelman et al.
# 2016 Current Biology to generate "initial" (draft) haplotypes for
# estimating introgressed haplotype frequencies.
#Optional arguments:
# prefix: Prefix for the cluster IDs, e.g. population examined in S'
# debug:  Output debugging and diagnostic information to STDERR
#         e.g. the nodes and edges of the graphs to diagnose cases
#              where edges are missing
# thresh: Minimum r^2 required to count an edge (default: 0.3)
#         This argument should be unnecessary, as the input should be
#         pre-filtered using --min-r2 to decrease processing time and
#         file size
BEGIN{
   FS="\t";
   OFS=FS;
   PROCINFO["sorted_in"]="@ind_num_asc";
   if (length(thresh) == 0) {
      thresh=0.3;
   };
}
#Set up a hash for the columns to avoid issues if columns are shuffled:
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
#Keep track of the chromosome for input warnings and as a component of
# the cluster IDs:
NR==2{
   chrom=$cols["CHR"];
}
#Compile the cluster map and edge adjacency matrix/map:
NR>1{
   if ($cols["CHR"] != chrom) {
      print "This input file contains more than one chromosome/scaffold ("chrom", "$cols["CHR"]")" > "/dev/stderr";
      print "We don't handle inter-scaffold haplotype clustering." > "/dev/stderr";
      exit 2;
   };
   if (length($cols["POS1"]) == 0) {
      print "POS1 missing for line "NR > "/dev/stderr";
   };
   if (length($cols["POS2"]) == 0) {
      print "POS2 missing for line "NR > "/dev/stderr";
   };
   if (length($cols["R^2"]) == 0) {
      print "R^2 missing for line "NR > "/dev/stderr";
   };
   if ($cols["R^2"] >= thresh) {
      if (($cols["POS1"] in nodes) && ($cols["POS2"] in nodes)) {
         if (length(nodes[$cols["POS1"]]) == 0) {
            print "Somehow both SNPs are in nodes, but one's cluster is blank..." > "/dev/stderr";
            print $0 > "/dev/stderr";
         };
         if (nodes[$cols["POS1"]] != nodes[$cols["POS2"]]) {
            #This is the main error case
            if (length(debug) > 0) {
               print "Edge potentially connecting clusters:" > "/dev/stderr";
               print chrom, $cols["POS1"], $cols["POS2"], $cols["R^2"], nodes[$cols["POS1"]], nodes[$cols["POS2"]] > "/dev/stderr";
            };
         } else {
            #Totally normal case, don't do anything as the nodes
            # are already in the same cluster, we just need to
            # set the edge.
         };
      } else if ($cols["POS1"] in nodes) {
         #The second SNP hasn't been seen before, so this is a normal
         # edge, though not the first in the cluster. Set the second
         # SNP to the same cluster:
         nodes[$cols["POS2"]]=nodes[$cols["POS1"]];
      } else if ($cols["POS2"] in nodes) {
         #The first SNP hasn't been seen before, so there's probably
         # a missing edge in this cluster. Set the first SNP to the
         #same cluster:
         nodes[$cols["POS1"]]=nodes[$cols["POS2"]];
      } else {
         #Neither SNP has been seen before, so this is a new cluster:
         nodes[$cols["POS1"]]=$cols["POS1"];
         nodes[$cols["POS2"]]=nodes[$cols["POS1"]];
      };
      #Add the undirected edge:
      #(vcftools output is necessarily in i < j order, so we only
      # need to represent one orientation)
      edges[$cols["POS1"],$cols["POS2"]]=1;
   };
}
#In post-processing, we'll rename the clusters and establish their
# completeness:
END{
   #Count missing edges within each cluster:
   for (i in nodes) {
      if (!(nodes[i] in missing)) {
         missing[nodes[i]]=0;
      };
      clustersize[nodes[i]]++;
      for (j in nodes) {
         if (i < j) {
            if (!((i,j) in edges) && nodes[i] == nodes[j]) {
               ++missing[nodes[i]];
            };
         };
      };
   };
   #Rename the clusters and output each SNP in order:
   print "CHROM", "POS", "CLUSTERID", "CLUSTERSIZE", "MISSINGEDGES";
   cluster_index=0;
   for (p in nodes) {
      if (!(nodes[p] in idmap)) {
         idmap[nodes[p]]=prefix"_"chrom"_"++cluster_index;
      };
      print chrom, p, idmap[nodes[p]], clustersize[nodes[p]], missing[nodes[p]];
   };
   #Diagnostic info:
   if (length(debug) > 0) {
      print "CHROM", "NODE/POS", "CLUSTERPOS" > "/dev/stderr";
      for (p in nodes) {
         print chrom, p, nodes[p] > "/dev/stderr";
      };
      print "CHROM", "NODE1/POS1", "NODE2/POS2", "CLUSTERPOS1", "CLUSTERPOS2", "HASEDGE" > "/dev/stderr";
      for (i in nodes) {
         for (j in nodes) {
            if (i < j) {
               print chrom, i, j, nodes[i], nodes[j], length(edges[i,j]) > "/dev/stderr";
            };
         };
      };
   };
}
