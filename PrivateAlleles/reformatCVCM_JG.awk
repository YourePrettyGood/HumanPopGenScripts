#!/bin/awk -f
#The filename parsing is custom for the joint genotyping CVCM files
#Required arguments:
# prefix:        Prefix of the filename prior to the dbSNP version
#                (I used "PIB92_SD_hs37d5_SNP_VQSR99.5_INDEL_VQSR99.0_dbSNP")
# suffix:        Suffix of the filename after the population label
#                (I used "_CVCM.variant_calling_summary_metrics")
# intermediate1: String between dbSNP version and subsampling ID
#                (I used "_VQSRpass_")
# intermediate2: String between subsampling ID and chromosome number
#                (I used "_chr")
# intermediate3: String between chromosome number and localization of site
#                (I used "_")
BEGIN{
   FS="\t";
   OFS=FS;
   ENF=25;
   filenum=0;
}
FNR==1{
   filenum++;
   n_path_parts=split(FILENAME, pathparts, "/");
   fn=pathparts[n_path_parts];
   sub(prefix, "", fn);
   sub(suffix, "", fn);
   sub(intermediate1, "/", fn);
   split(fn, a, "/");
   dbsnpversion=a[1];
   fn=a[2];
   sub(intermediate2, "/", fn);
   split(fn, a, "/");
   subsampling=a[1];
   fn=a[2];
   sub(intermediate3, "/", fn);
   split(fn, a, "/");
   chrom=a[1];
   fn=a[2];
   if (fn ~ /^subset/) {
      sub("subset", "", fn);
      region=fn;
      sites="all";
   } else {
      sub("private", "", fn);
      region=fn;
      sites="private";
   };
   skip=1;
}
#Feed through the statistic lines, prepending some key columns:
skip==0&&NF==ENF{
   print chrom, dbsnpversion, subsampling, region, sites, $0;
}
#SAMPLE_ALIAS is the first column header in the detail_metrics files
#TOTAL_SNPS is the first column header in the summary_metrics files
/^SAMPLE_ALIAS/||/^TOTAL_SNPS/{
   if (filenum==1) {
      print "Chromosome", "dbSNPversion", "Subsamples", "Region", "Sites", $0;
   };
   skip=0;
   ENF=NF;
}
