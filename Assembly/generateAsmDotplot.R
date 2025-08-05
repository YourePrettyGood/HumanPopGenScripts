#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, repos="https://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("tidyverse")

#Read in arguments:
paf_fn <- options[1]
ctg_order_fn <- options[2]
ref_mask_fn <- options[3]
query_id <- options[4]
ref_id <- options[5]
plot_fn <- options[6]

#Load in the inferred order of query assembly contigs and corresponding ref chromosomes:
ctg_order <- read_tsv(ctg_order_fn,
                      col_names=c("Query", "Ref", "OrderIndex", "gpOrientation", "Orientation"),
                      col_types='ccicc')
#Load in the PAF of alignments between query assembly contigs and ref chromosomes:
#The col_names here only handles the default minimap2 output tags and will fail
# if there are more or less than 5 tags on a given line.
paf <- read_tsv(paf_fn,
                col_names=c("Query", "QueryLen", "QueryStartBED", "QueryEnd", "AlnOrientation",
                            "Ref", "RefLen", "RefStartBED", "RefEnd", "Match", "AlnLen",
                            "MQ", "tagA", "tagB", "tagC", "tagD", "tagE"),
                col_types='cnnnccnnnnniccccc')

#Extract the lengths of the query contigs and ref chromosomes:
query_lens <- paf %>%
   distinct(Query, QueryLen)
ref_lens <- paf %>%
   distinct(Ref, RefLen) %>%
   filter(Ref != "*")

#Determine the offsets of each contig and chromosome along the plot axes:
#This is done because the plot axes are along a single coordinate, so query
# and ref alignment positions need to be offset along this coordinate to
# visually correspond to their query contig and ref chromosome.
ctg_offsets <- ctg_order %>%
   right_join(query_lens, by="Query") %>%
   transmute(Query=Query,
             QueryOffset=lag(cumsum(QueryLen), default=0))
chr_offsets <- ctg_order %>%
   distinct(Ref) %>%
   right_join(ref_lens, by="Ref") %>%
   transmute(Ref=Ref,
             RefOffset=lag(cumsum(RefLen), default=0))

#Read in the BED of masked regions (i.e. centromeres, telomeres, PARs):
#These regions get coloured transparent brown to indicate where query
# contigs might have expected breaks (and where some query contigs might
# unexpectedly traverse).
mask_df <- read_tsv(ref_mask_fn,
                    col_names=c("Ref", "MaskStartBED", "MaskEndBED"),
                    col_types='cnn') %>%
   inner_join(chr_offsets, by="Ref") %>%
   transmute(Ref=Ref,
             MaskStart=MaskStartBED+RefOffset,
             MaskEnd=MaskEndBED+RefOffset-1,
             QueryStart=0,
             QueryEnd=sum(query_lens$QueryLen))

#Reformat the PAF so that alignments are now in plot coordinate space
# by using the offsets we calculated earlier:
#Make sure to adjust coordinates based on alignment orientation.
aln_df <- paf %>%
   left_join(ctg_offsets, by="Query") %>%
   left_join(chr_offsets, by="Ref") %>%
   transmute(Query=Query,
             QueryStart=QueryOffset+QueryStartBED,
             QueryEnd=QueryOffset+QueryEnd-1,
             AlnOrientation=AlnOrientation,
             Ref=Ref,
             RefStart=case_when(AlnOrientation == "-" ~ RefEnd+RefOffset-1,
                                TRUE ~ RefStartBED+RefOffset),
             RefEnd=case_when(AlnOrientation == "-" ~ RefStartBED+RefOffset,
                              TRUE ~ RefEnd+RefOffset-1),
             QueryOffset=QueryOffset,
             RefOffset=RefOffset)

#Now generate the actual dotplot in the MUMmer/mummerplot style, with
# points at each end of each alignment and a segment between the ends,
# coloured red if + orientation and blue if - orientation:
#Grey dashed lines indicate the boundaries between contigs/chromosomes.
#Contigs and chromosomes are labeled at their starting coordinates.
dotplot <- aln_df %>%
   filter(AlnOrientation != "*") %>%
   mutate(AlnOrientation=factor(AlnOrientation, levels=c("+", "-"))) %>%
   ggplot() +
      geom_rect(data=mask_df,
                aes(xmin=MaskStart,
                    xmax=MaskEnd,
                    ymin=QueryStart,
                    ymax=QueryEnd),
                fill="brown",
                alpha=0.3,
                inherit.aes=FALSE) +
      geom_hline(aes(yintercept=QueryOffset),
                 col="grey50",
                 alpha=0.3,
                 linetype=3) +
      geom_vline(aes(xintercept=RefOffset),
                 col="grey50",
                 alpha=0.7,
                 linetype=3) +
      geom_point(aes(x=RefStart,
                     y=QueryStart,
                     colour=AlnOrientation)) +
      geom_point(aes(x=RefEnd,
                     y=QueryEnd,
                     colour=AlnOrientation)) +
      geom_segment(aes(x=RefStart,
                       xend=RefEnd,
                       y=QueryStart,
                       yend=QueryEnd,
                       colour=AlnOrientation)) +
      theme_classic() +
      scale_colour_manual(values=c("red", "blue")) +
      scale_x_continuous(name=ref_id,
                         breaks=chr_offsets$RefOffset,
                         labels=chr_offsets$Ref,
                         expand=c(0,0)) +
      scale_y_continuous(name=query_id,
                         breaks=ctg_offsets$QueryOffset,
                         labels=ctg_offsets$Query,
                         expand=c(0,0)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#Save the plot as a 24.0x24.0 cm PDF at 500 DPI:
ggsave(plot_fn,
       plot=dotplot,
       width=24.0,
       height=24.0,
       units="cm",
       dpi=500)
