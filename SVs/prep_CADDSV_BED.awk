#!/bin/awk -f
#This script takes the TSV output by select_SVs_wAFs.awk and
# converts it to a BED file needed as input for CADD-SV.
#The CADD-SV BED format involves the standard first three
# BED columns (chromosome, 0-based inclusive start, 0-based
# exclusive end aka 1-based inclusive end),
# as well as an SVTYPE column and an SV sequence column.
#The SVTYPE column should be DEL, DUP, INS, or INV.
#The SV sequence column should be empty if DEL, the
# inserted sequence if INS, and I have no idea for
# DUP or INV.
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!("CHROM" in cols)) {
      print "Missing CHROM column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (!("POS" in cols)) {
      print "Missing POS column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 3;
   };
   if (!("SVTYPE" in cols)) {
      print "Missing SVTYPE column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 4;
   };
   if (!("SVLEN" in cols)) {
      print "Missing SVLEN column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 5;
   };
   if (!("REF" in cols)) {
      print "Missing REF column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 6;
   };
   if (!("ALT" in cols)) {
      print "Missing ALT column from input, cannot proceed. Quitting." > "/dev/stderr";
      exit 7;
   };
}
NR>1{
   if ($cols["SVTYPE"] == "INS") {
      bedstart=$cols["POS"]-1;
      bedend=$cols["POS"];
      refseq=$cols["REF"];
      altseq=$cols["ALT"];
      bedseq=substr(altseq, length(refseq)+1);
      if (length(bedseq) != $cols["SVLEN"]) {
         if (length(debug) > 0) {
            print "Somehow length of insertion sequence (altlen "length(altseq)" - reflen "length(refseq)" = "length(bedseq)") is not equal to SVLEN ("$cols["SVLEN"]") for "$cols["SVTYPE"]" at "$cols["CHROM"]":"$cols["POS"] > "/dev/stderr";
         };
         next;
      };
   } else if ($cols["SVTYPE"] == "DEL") {
      #Don't subtract one for start, since DEL has one base padding on left
      # in VCF:
      bedstart=$cols["POS"];
      bedend=bedstart+$cols["SVLEN"];
      bedseq="";
   } else {
      if (length(debug) > 0) {
         print "Unsupported SVTYPE "$cols["SVTYPE"]" for the time being." > "/dev/stderr";
      };
      next;
   };
   if (length(svid) > 0) {
      print $cols["CHROM"], bedstart, bedend, $cols["SVTYPE"], bedseq, $cols["CHROM"]"-"$cols["POS"]"-"$cols["SVTYPE"]"-"$cols["SVLEN"];
   } else {
      print $cols["CHROM"], bedstart, bedend, $cols["SVTYPE"], bedseq;
   };
}
