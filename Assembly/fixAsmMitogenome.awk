#!/bin/awk -f
#This script takes a base-level alignment PAF from minigraph or minimap2 of
# a query assembly against a reference genome as well as the query assembly
# FASTA, and identifies candidates for a full-length mitogenome contig,
# then trims and revcomps the best candidate, and drops the rest.
#The output is a modified version of the query assembly FASTA with the
# top mitogenome contig candidate trimmed and renamed (with an "m" suffix
# instead of hifiasm's typical "l" or "c") and all other mitogenome
# candidate contigs dropped.
#We use a minimum coverage threshold for the candidate list to avoid dropping
# contigs with NUMTs, which tend to be a few kbp rather than the full
# (human) 16.5 kbp. The default mincov is 0.7, meaning at least 70% coverage
# of the reference mitogenome must be covered by the alignment (including
# gaps). Technically, this isn't precisely "coverage", it's just the ratio
# of the number of columns (including gaps) in the alignment to the length
# of the reference mitogenome, but it should be pretty close.
#I don't handle extreme cases here, since the expectation is that your input
# is hifiasm contigs from e.g. humans with a human reference genome, so the
# divergence of your candidate mitogenome contigs to the reference should
# be sufficiently low to work with this heuristic.
#The top candidate is identified as the query contig that contains an
# alignment with gapped coverage the closest to 1.0 of all candidates.
#This is achieved by calculating the absolute value of the reference
# coverage deviation from 1.0 and taking the minimum.
#Options:
# mincov: Minimum reference mitogenome coverage by alignment for a contig
#         to qualify as a candidate (default: 0.7)
# mtid:   ID of the mitogenome contig in the reference genome
#         (default: chrM)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(mincov) == 0) {
      mincov=0.7;
      print "No minimum mitogenome coverage threshold provided, using 0.7." > "/dev/stderr";
   };
   if (mincov < 0.0 || mincov > 1.0) {
      print "Invalid mincov provided, must be between 0 and 1. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(mtid) == 0) {
      mtid="chrM";
   };
   n_candidates=0;
   #Flag to extract the current contig:
   extract=0;
   #Flag to drop the current contig:
   drop=0;
   #Keep track of which file we're on:
   filenum=0;
   #Hash for transliterating complementary bases:
   #Start with the critical set of complementing definitions:
   #We skip S, W, and N, since they're just idempotent.
   critbase["A"]="T";
   critbase["C"]="G";
   critbase["R"]="Y";
   critbase["K"]="M";
   critbase["B"]="V";
   critbase["D"]="H";
   #Account for case and symmetry:
   for (base in critbase) {
      compbase[base]=critbase[base];
      compbase[tolower(base)]=tolower(critbase[base]);
      compbase[critbase[base]]=base;
      compbase[tolower(critbase[base])]=tolower(base);
   };
}
function revcomp(seq, complhash) {
   rcseq="";
   n=split(seq, bases, "");
   for (i=n; i>0; i--) {
      if (bases[i] in complhash) {
         rcseq=rcseq complhash[bases[i]];
      } else {
         rcseq=rcseq bases[i];
      };
   };
   return rcseq;
}
#Keep track of which file we're on:
FNR==1{
   filenum+=1;
}
#First file is the PAF of base-level alignment to ref genome.
#We first identify candidate near-complete or complete alignments to the
# mtDNA contig in the ref genome and store coordinates and orientation:
filenum==1&&$6==mtid{
   if ($11 >= ($7 * mincov)) {
      n_candidates+=1;
      #Store the query contig name of the candidate:
      candidate[n_candidates]=$1;
      #And keep them in a hash so that the rest can be dropped:
      candidates[$1]=n_candidates;
      #The unsigned deviation of reference coverage of the alignment from 1.0
      #(to sort later for prioritizing candidates):
      signedcovdeviation=1.0-($11/$7);
      alncovdeviation[n_candidates]=signedcovdeviation < 0.0 ? -1 * signedcovdeviation : signedcovdeviation;
      #The 1-based start of the alignment on the query contig (to use with substr):
      qstart[n_candidates]=$3+1;
      #The length of the alignment in query contig space (to use with substr):
      qalnlen[n_candidates]=$4-$3;
      #The orientation of the alignment (i.e. should we revcomp?):
      orient[n_candidates]=$5;
   };
}
#Second file is the FASTA to parse to retrieve the mtDNA contig.
#On the first line, determine the index of the top candidate:
filenum==2&&FNR==1{
   PROCINFO["sorted_in"]="@val_num_asc";
   for (qidx in alncovdeviation) {
      topq=qidx;
      print "Top mtDNA candidate contig is "candidate[topq] > "/dev/stderr";
      break;
   };
}
#Now check the header lines to find the appropriate contig and output it:
filenum==2&&/^>/{
   #Trim, revcomp as needed, and output the previous contig if requested:
   if (extract == 1) {
      #Modify the hifiasm contig ID with an "m" suffix instead of "l" or "c" to indicate mtDNA:
      if (ctgid ~ /[lc]$/) {
         print ">"substr(ctgid, 1, length(ctgid)-1)"m";
      } else {
         print ">"ctgid"m";
      };
      if (orient[topq] == "-") {
         print revcomp(substr(seq, qstart[topq], qalnlen[topq]), compbase);
      } else {
         print substr(seq, qstart[topq], qalnlen[topq]);
      };
   } else if (drop == 0 && FNR>1) {
      print ">"ctgid;
      print seq;
   };
   #Reset the extract and drop states:
   extract=0;
   drop=0;
   #Also reset the stored sequence:
   seq="";
   #Get rid of the > prefix so that we can test contig ID equality:
   ctgid=substr($1, 2);
   #Determine if we need to do extraction or dropping of the current contig:
   if (ctgid == candidate[topq]) {
      extract=1;
   } else if (ctgid in candidates) {
      drop=1;
   };
}
#Store the sequence of the best mtDNA candidate contig:
filenum==2&&!/^>/{
   if (extract == 1 || drop == 0) {
      seq=seq $0;
   };
}
END{
   #Catch the case where the contig to extract or not drop was the last one in the FASTA:
   if (extract == 1) {
      #Modify the hifiasm contig ID with an "m" suffix instead of "l" or "c" to indicate mtDNA:
      if (ctgid ~ /[lc]$/) {
         print ">"substr(ctgid, 1, length(ctgid)-1)"m";
      } else {
         print ">"ctgid"m";
      };
      if (orient[topq] == "-") {
         print revcomp(substr(seq, qstart[topq], qalnlen[topq]), compbase);
      } else {
         print substr(seq, qstart[topq], qalnlen[topq]);
      };
   } else if (drop == 0 && FNR>1) {
      print ">"ctgid;
      print seq;
   };
}
