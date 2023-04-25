#!/bin/awk -f
#This script takes the output of bcftools query -f '%CHROM:%POS\t%REF\t%ALT\n'
# on the MAF-pruned VCF of the OAS region as well as the CADD *.rankRawScore.gz
# file for the appropriate chromosome to add a CADD score column to each
# site.
#Note that these CADD scores must be deduped by taking the maximum
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum++;
}
#The first file is the output of bcftools query -f '%CHROM:%POS\t%REF\t%ALT\n'
# on the VCF of sites in the region, so we store REF and ALT:
filenum==1{
   site[$1]=$2 SUBSEP $3;
}
#The second file is the decompressed CADD *.rankRawScore.gz file from WGSA:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Output a header:
   print "CHROM:POS", "REF", "ALT", "CADD";
}
filenum==2&&FNR>1{
   sitekey=$1":"$2;
   if (sitekey in site) {
      split(site[sitekey], refalt, SUBSEP);
      n_alt=split(refalt[2], alts, ",");
      for (j=1; j<=n_alt; j++) {
         if ($cols["Ref"] == refalt[1] && $cols["Alt"] == alts[j]) {
            print sitekey, $cols["Ref"], $cols["Alt"], $cols["PHRED"];
         };
      };
   };
}
