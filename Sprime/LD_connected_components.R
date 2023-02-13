#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, repos="https://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("igraph")
check_package("tidyverse")

#Read in arguments:
input <- options[1]
output <- options[2]
chrom <- options[3]

#Read in the edges and convert to a data.frame compatible with igraph:
LD <- read_tsv(input) %>%
   transmute(POS_A=POS_A, POS_B=POS_B, Rsq=Rsquared)
#Create a graph from the edges in the data.frame and all inferred nodes:
LD_graph <- graph_from_data_frame(LD, directed=FALSE)
#Determine the connected components in the graph:
#Note: mode is ignored, since the graph is undirected.
LD_CCs <- components(LD_graph, mode="strong")
#Quick diagnostic of the number of connected components found
LD_CCs$no
#Reformat the components into lists of tag SNPs and their haplotype
# IDs:
LD_clusters <- data.frame(CHROM=chrom,
                          POS=as.integer(names(LD_CCs$membership)),
                          ClusterID=LD_CCs$membership) %>%
   arrange(CHROM, ClusterID)
#Write the inferred core haplotypes to a file:
write_tsv(LD_clusters,
          file=output,
          quote="none")
