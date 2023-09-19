#!/bin/awk -f
#This script takes the output of bctools query -H -f '%CHROM:%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n'
# on a phased VCF and outputs one haplotype per column in TSV format.
#Ancestral allele index is included as the second column if the ancestral allele is input.
#The output can be used easily in R to plot haplotypes as a heatmap.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(ploidy) == 0) {
      ploidy=2;
      print "Defaulting to ploidy of 2" > "/dev/stderr";
   };
   #Set default for flag for ancestral allele in input:
   has_aa=0;
}
NR==1{
   header="";
   for (i=1; i<=NF; i++) {
      sub("#[ ]?", "", $i);
      gsub("[[0-9]+]", "", $i);
      if ($i ~ ":GT") {
         gsub(":GT", "", $i);
         gtcols[i]=1;
         for (j=1; j<=ploidy; j++) {
            header=header OFS $i"_"j;
         };
      } else {
         cols[$i]=i;
      };
   };
   if (("AA" in cols) && ("REF" in cols) && ("ALT" in cols)) {
      print "Found ancestral allele in column names, outputting ancestral allele index column" > "/dev/stderr";
      has_aa=1;
      header=OFS "Ancestral" header;
   };
   if ("CHROM:POS" in cols) {
      header="CHROM:POS" header;
   } else {
      print "Input is missing CHROM:POS column, cannot proceed." > "/dev/stderr";
      exit 2;
   };
   print header;
}
NR>1{
   line=$cols["CHROM:POS"];
   if (has_aa > 0) {
      anc_allele="NA";
      n_alt=split($cols["ALT"], alleles, ",");
      alleles[0]=$cols["REF"];
      for (i=0; i<=n_alt; i++) {
#         print $cols["CHROM:POS"], $cols["REF"], $cols["ALT"], $cols["AA"], alleles[i] > "/dev/stderr";
         if ($cols["AA"] ~ /^[ACGTacgt]$/ && $cols["AA"] == alleles[i]) {
            anc_allele=i;
         };
      };
      line=line OFS anc_allele;
   };
   for (i=1; i<=NF; i++) {
      if (i in gtcols) {
         n=split($i, alleles, "[/|]");
         for (j=1; j<=ploidy; j++) {
            sub("[.]", "NA", alleles[j]);
            line=line OFS alleles[j];
         };
      };
   };
   print line;
}
