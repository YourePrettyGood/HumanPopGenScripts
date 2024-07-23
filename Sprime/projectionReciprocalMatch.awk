#!/bin/awk -f
#This script takes BEDs produced by addOriginToPhasedProjections.awk
# as well as the output of bedtools intersect -wo between
# pairs of these BEDs and summarizes the values of
# |X_i \intersect X_j| and |X_i| for individuals
# i and j (where i comes from file a and j comes from file b
# of the intersect call).
#The value of |X_i \intersect X_j| is simply the count of
# overlapping projections between the two individuals.
#We can use these values along with the total numbers of
# projections for each individual (|X_i|) to calculate the
# conditional probability of archaic tract sharing.
#Thus, the output contains estimates of Pr(i|j) and Pr(j|i).
#It also contains all necessary constituents to recalculate
# Pr(i|j) and Pr(j|i) in another language in case extra precision
# is desired.
#Mandatory arguments:
#  origin:     Archaic origin of the input tracts (e.g. Neanderthal,
#              Denisovan, Ambiguous)
#Optional arguments:
#  header:     Whether or not to print a header line for the output
#              (default: no header)
#  minoverlap: Minimum length of overlap (in bp)  between projections to count
#              (default: 1)
BEGIN{
   FS="\t";
   OFS=FS;
   #Default minimum overlap length to count is 1 bp:
   if (length(minoverlap) == 0) {
      minoverlap=1;
   };
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum+=1;
};
#First file is the concatenated output of addOriginToPhasedProjections.awk:
#We get the values of \hat{Pr(i)} from this file.
#Skip any header lines and identify i to add to its count:
filenum==1&&!/^#/{
   #Typical layout of addOriginToPhasedProjections.awk output:
   # 1) Scaffold/Chromosome
   # 2) BED-style haplotype start position (i.e. 0-based)
   # 3) BED-style haplotype end position (i.e. 1-based)
   # 4) Sprime tract ID (i.e. [Population]_[Chromosome]_[Segment ID])
   # 5) Haplotype ID (i.e. [Individual ID]_[Haplotype #])
   # 6) Haplotype state (ignored, has issues/bugs)
   # 7) # Sprime-matching sites on haplotype
   # 8) # Sprime-mismatching sites on haplotype
   # 9) Neanderthal match rate
   # 10) Denisovan match rate
   # 11) Altai Neanderthal match rate
   # 12) Vindija Neanderthal match rate
   # 13) Chagyrskaya Neanderthal match rate
   #So we only care about column 5 from this file.
   #Check if the input seems to have the correct number of columns:
   if (NF != 13) {
      print "Is the first file ("FILENAME") actually the output of addOriginToPhasedProjections.awk? It doesn't look like it, as there are "NF" columns, not 13." > "/dev/stderr";
      exit 2;
   };
   if (length($5) == 0) {
      print "Haplotype ID is empty for line "NR > "/dev/stderr";
   };
   hapid=$5;
   jointprob[hapid,hapid]+=1;
}
#Second file is the concatenated output of bedtools intersect -wo
# between populations:
#We get the values of \hat{Pr(i&j)} from this file.
#Skip any header lines and identify i and j to add to their count:
filenum==2&&!/^#/{
   #Check if the input seems to have the correct number of columns:
   if (NF != 27 && NF != 28) {
      print "Is the first file ("FILENAME") actually the output of bedtools intersect -wo? It doesn't look like it, as there are "NF" columns, not 27 or 28." > "/dev/stderr";
      exit 3;
   };
   #Determine the offset between the A and B columns of the bedtools
   # intersect -wo output:
   if (NF % 2 == 1) {
      #When -b is an individual BED, NF will be 27, since there will be
      # 13 columns per input BED followed by a single column for the
      # length of overlap, so 2*13+1=27.
      #Therefore, the offset will be (27-1)/2=13
      offset=(NF-1)/2;
   } else {
      #When -b is a list of BEDs, NF will be 28, since there will be
      # 13 columns for the query BED, a column indicating which
      # target BED the overlap is from, then 13 columns for the
      # target BED, and a final column for the length of overlap,
      # so 13+1+13+1=28.
      #Therefore, the offset will be 28/2=14
      offset=NF/2;
   };
   #Extract the length of the overlap:
   overlap=$NF;
   if (overlap >= minoverlap) {
      #Typical layout of addOriginToPhasedProjections.awk output:
      # 1) Scaffold/Chromosome
      # 2) BED-style haplotype start position (i.e. 0-based)
      # 3) BED-style haplotype end position (i.e. 1-based)
      # 4) Sprime tract ID (i.e. [Population]_[Chromosome]_[Segment ID])
      # 5) Haplotype ID (i.e. [Individual ID]_[Haplotype #])
      # 6) Haplotype state (ignored, has issues/bugs)
      # 7) # Sprime-matching sites on haplotype
      # 8) # Sprime-mismatching sites on haplotype
      # 9) Neanderthal match rate
      # 10) Denisovan match rate
      # 11) Altai Neanderthal match rate
      # 12) Vindija Neanderthal match rate
      # 13) Chagyrskaya Neanderthal match rate
      #Thus, offset should be 13 unless the code has changed.
      qscaf=$1;
      qstart=$2;
      qend=$3;
      qtract=$4;
      qhap=$5;
      tscaf=$(offset+1);
      tstart=$(offset+2);
      tend=$(offset+3);
      ttract=$(offset+4);
      thap=$(offset+5);
      #Skip non-overlap records in case -wao was accidentally used
      # instead of -wo:
      if (tstart < 0) {
         next;
      };
      #Also skip any records that were accidentally included where j==i:
      #These would mangle the counts we got from the original BEDs for Pr(i).
      if (qhap == thap) {
         next;
      };
      #We want to add this tract pair to \hat{Pr(i&j)}:
      #Note: i and j are haplotype IDs
      #Also, this isn't technically \hat{Pr(i&j)}, it's unnormalized,
      # but when divided by the unnormalized analog of \hat{Pr(i)},
      # we get \hat{Pr(i|j)} as specified in the supplement of
      # Vernot et al. 2016 Science.
      jointprob[qhap,thap]+=1;
   };
}
END{
   #Set the deterministic array iteration order:
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Output the counts:
   if (length(header) > 0 && header > 0) {
      print "origin", "i", "j", "Pr(i,j)", "Pr(i)", "Pr(j)", "Pr(i|j)", "Pr(j|i)";
   };
   for (pair in jointprob) {
      split(pair, ij, SUBSEP);
      if (ij[1] == ij[2]) {
         continue;
      };
      if (!((ij[1],ij[1]) in jointprob)) {
         print "Pr("ij[1]") is missing, check the first input file." > "/dev/stderr";
         exit 4;
      };
      if (!((ij[2],ij[2]) in jointprob)) {
         print "Pr("ij[2]") is missing, check the first input file." > "/dev/stderr";
         exit 5;
      };
      pr_iandj=jointprob[pair];
      pr_i=jointprob[ij[1],ij[1]];
      pr_j=jointprob[ij[2],ij[2]];
      pr_igivenj=pr_iandj/pr_i;
      pr_jgiveni=pr_iandj/pr_j;
      print origin, ij[1], ij[2], pr_iandj, pr_i, pr_j, pr_igivenj, pr_jgiveni;
   };
}
