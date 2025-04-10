#!/bin/awk -f
#Initially I'm only writing to outupt SNPs
#Deletions will be treated as missing SNP genotypes
#Insertions will be treated as a SNP at the flanking base and the
# inserted sequence will be ignored
#It should be possible to reconstruct indels from the EMF alignment,
# but I'm just writing this to get some quick results.
#
#Inputs:
# First file is the reference genome FAI
# Second file is the uncompressed EMF
# ref: ID of genome to take as reference
# query: ID of genome to determine genotype for
# progress: Output a progress indicator every x alignment blocks (default: 100)
# haploid: Force the genotype to be output as haploid
# qual: Value of QUAL to put in the output (default: 30)
# filter: Value of FILTER to put in the output (default: .)
#Define functions for reverse-complementing alleles:
function complement(x) {
   #We don't cover IUPAC degenerate bases for now.
   X=toupper(x);
   if (X == "A") {
      xbar="T";
   } else if (X == "C") {
      xbar="G";
   } else if (X == "G") {
      xbar="C";
   } else if (X == "T") {
      xbar="A";
   } else if (X == "N") {
      xbar="N";
   } else {
      xbar=".";
   };
   if (x != X) {
      xbar=tolower(xbar);
   };
   return xbar;
}
function revcomp(allele) {
   rc_allele="";
   allele_len=split(allele, seq, "");
   for (i=allele_len; i>=1; i--) {
      rc_allele=rc_allele""complement(seq[i]);
   };
   return rc_allele;
}
BEGIN{
   version="1.0";
   OFS="\t";
   cmdline="-v ref="ref" -v query="query;
   if (length(haploid) > 0) {
      cmdline=cmdline" -v haploid="haploid;
   };
   if (length(progress) == 0) {
      progress=100;
   } else {
      cmdline=cmdline" -v progress="progress;
   };
   if (length(qual) == 0) {
      qual=30;
   } else {
      cmdline=cmdline" -v qual="qual;
   };
   print "Outputting VCF records with QUAL forced to "qual > "/dev/stderr";
   if (length(filter) == 0) {
      filter=".";
   } else {
      cmdline=cmdline" -v filter="filter;
   };
   cmdline=cmdline" "ARGV[1]" "ARGV[2];
   print "Outputting VCF records with FILTER forced to "filter > "/dev/stderr";
   readin=0;
   seqid=0;
   blockid=0;
   skipBlock=0;
   #For keeping track of the input file:
   filenum=0;
   #Beginning of the VCF header:
   print "##fileformat=VCFv4.2";
   print "##source=EMFtoVCFheader.awk version "version;
   print "##EMFtoVCFheaderLine=<ID=EMFtoVCFheader,CommandLine=\"EMFtoVCFheader.awk "cmdline"\">";
   #Major caveat to this: One must pass the absolute path to the FAI
   # in order to produce a valid URI here:
   refpath=ARGV[1];
   sub(/.fai$/, "", refpath);
   print "##reference=file:///"refpath;
}
FNR==1{
   filenum++;
}
filenum==1{
   #Print the reference scaffold lines of the VCF header:
   print "##contig=<ID="$1",length="$2">";
}
filenum>1&&FNR==1{
   #Print INFO and FORMAT metadata lines of the VCF header:
   print "##INFO=<ID=EMFBLOCK,Number=1,Type=Integer,Description=\"0-based alignment block index from the EMF passed to EMFtoVCFheader.awk\",Source=\"EMFtoVCFheader.awk\",Version=\""version"\">";
   print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
   #Print the #CHROM line of the VCF header:
   print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", query;
}
filenum>1&&/^SEQ/{
   #Skip any blocks that aren't 1:1 alignments, as we can't represent those
   # in VCF format:
   if ($2 == ref && ref in colids) {
      print "Refdup: Skipping non-1-to-1 alignment block "blockid > "/dev/stderr";
      skipBlock=1;
   };
   if ($2 == query && query in colids) {
      print "Querydup: Skipping non-1-to-1 alignment block "blockid > "/dev/stderr";
      skipBlock=1;
   };
   #Collect info about the sequences in this block, including a mapping of
   # node positions in the alignment column:
   #Note: One failure point here is that Ensembl's EMFs use the same genome
   # ID for each of the ancestral reconstructions, so the bug here is that
   # the first ancestral reconstruction in the alignment is always output.
   colids[$2]=seqid;
   ids[seqid]=$2;
   chroms[seqid]=$3;
   starts[seqid]=$4;
   ends[seqid]=$5;
   strands[seqid]=$6 > 0 ? "+" : "-";
   seqid++;
}
filenum>1&&/^DATA/{
   readin=1;
   alnid=0;
   refskip=0;
   #We might want to do something to deal with ambiguous ancestral state
   # IDs, as ids will contain "ancestor_sequences" multiple times
   #For now we don't do anything.
   #Make sure the DATA line doesn't get interpreted by later rules:
   next;
}
filenum>1&&!/^\/\//&&readin==1{
   #Skip the block when requested:
   if (skipBlock) {
      next;
   };
   #Split the current alignment column into an array, one element per node:
   n_nodes=split($0, states, "");
   #Check to make sure the EMF isn't malformatted:
   if (n_nodes != seqid) {
      print "Number of nodes in alignment column ("n_nodes") does not match the number of sequences in the alignment ("seqid") for block "blockid" with alignment column:" > "/dev/stderr";
      print $0 > "/dev/stderr";
      exit 2;
   };
   #Now check the ref and query states:
   mask=0;
   refstate=states[colids[ref]+1];
   querystate=states[colids[query]+1];
   if (toupper(refstate) ~ /[ACGT]/) {
      if (refstate != toupper(refstate)) {
         mask=1;
      };
      alleles[0]=toupper(refstate);
      #non-gap:
      if (toupper(querystate) ~ /[ACGT]/) {
         #Check for match or mismatch:
         if (toupper(querystate) != toupper(refstate)) {
            alleles[1]=toupper(querystate);
            gt="1/1";
         } else {
            alleles[1]=".";
            gt="0/0";
         };
      #Ref is defined, but query is not, so for now we treat this as missing
      # genotype in query, and no alt:
      } else {
         alleles[1]=".";
         gt="./.";
      };
      #Here's where we handle masking:
      if (mask) {
         alleles[1]=".";
         gt="./.";
      };
   #Ref is defined, but is non-gap and non-ACGT, so we treat it as missing
   # genotype in query, and no alt:
   } else if (refstate ~ /[^ACGTacgt~-]/) {
      alleles[0]=toupper(refstate);
      alleles[1]=".";
      gt="./.";
   #Ref is not defined, so for now we treat this as a skip, as it's a likely
   # insertion:
   } else {
      refskip++;
      alnid++;
      next;
   };
   #Adjust the alleles to account for ref strand:
   if (strands[colids[ref]] == "-") {
      alleles[0]=revcomp(alleles[0]);
      alleles[1]=revcomp(alleles[1]);
   };
   #Adjust the genotype for ploidy if requested:
   if (length(haploid) > 0) {
      split(gt, gtarr, "/");
      gt=gtarr[1];
   };
   #Now we calculate the position, and output the VCF line:
   if (strands[colids[ref]] == "-") {
      pos=ends[colids[ref]]-alnid+refskip;
#      print "row="alnid", refskip="refskip", boundary="ends[colids[ref]]", pos="pos > "/dev/stderr";
   } else {
      pos=starts[colids[ref]]+alnid-refskip;
#      print "row="alnid", refskip="refskip", boundary="starts[colids[ref]]", pos="pos > "/dev/stderr";
   };
   #Columns are CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO[, FORMAT, SAMPLE1[, SAMPLE...]]
   print chroms[colids[ref]], pos, ".", alleles[0], alleles[1], qual, filter, "EMFBLOCK="blockid, "GT", gt;
   #Keep track of the alignment column so we know how much to offset:
   alnid++;
}
filenum>1&&/^\/\//{
   readin=0;
   seqid=0;
   delete colids;
   delete ids;
   delete chroms;
   delete starts;
   delete ends;
   delete strands;
   skipBlock=0;
   blockid++;
   if (blockid % progress == 0) {
      print "Processed "blockid" blocks" > "/dev/stderr";
   };
}
