#!/bin/awk -f
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(threshcol) == 0) {
      print "No thresholding column name (threshcol) provided. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(thresh) == 0) {
      thresh=1;
      print "No threshold (thresh) value provided, using "threshcol" < 1 by default." > "/dev/stderr";
   };
   if (length(minlen) == 0) {
      minlen=0;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   if (!("CHROM" in cols)) {
      print "Missing CHROM column from input. Quitting." > "/dev/stderr";
      exit 3;
   };
   if (!("POS" in cols)) {
      print "Missing POS column from input. Quitting." > "/dev/stderr";
      exit 4;
   };
   if (!("ID" in cols)) {
      print "Missing ID column from input. Quitting." > "/dev/stderr";
      exit 5;
   };
   if (!(threshcol in cols)) {
      print "Missing "threshcol" column from input. Quitting." > "/dev/stderr";
      exit 6;
   };
   #Print the header:
   print "CHROM", "POS", "SVTYPE", "SVLEN", threshcol, "ID";
}
NR>1{
   #Only consider variants passing the threshold:
   if ($cols[threshcol] != "NA" && $cols[threshcol] != "." && $cols[threshcol] < thresh) {
      #Splitting ID assuming minigraph-cactus/PanGenie variant IDs:
      n=split($cols["ID"], idparts, "-");
      svtype=idparts[3];
      svlen=idparts[5];
      if (svlen >= minlen) {
         print $cols["CHROM"], $cols["POS"], svtype, svlen, $cols[threshcol], $cols["ID"];
      };
   };
}
