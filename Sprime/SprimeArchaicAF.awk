#!/bin/awk -f
#This script takes the .score file produced by S' and genotypes from a
# VCF processed by bcftools query -H -f '%CHROM:%POS[\t%GT]\n' and outputs
# the frequency of the putatively archaic allele identified by S' in the
# VCF.
#This is a preprocessing step necessary for identifying adaptively
# introgressed tracts, as we can take these per-site archaic allele
# frequencies and take the median AAF across sites in a tract as the
# tract frequency.
#Options:
# pop:    (required) The name of the populations the genotypes belong to
# header: (optional) A flag indicating whether or not to print the header
#                    line (default: don't print header)
#The "header" flag facilitates combining the results across chromosomes and
# populations.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   if (length(pop) == 0) {
      print "Missing pop variable, please set it." > "/dev/stderr";
      exit 2;
   };
   if (length(header) > 0) {
      print "CHROM", "POS", "TractID", "ArchaicAlleleCount", "TotalAlleleCount", "ArchaicAlleleFrequency", "Population";
   };
}
#Simple way to enable detection of more than 2 input files:
FNR==1{
   filenum++;
}
#First file is the S' score file, so store the column names for use later:
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      sprimecols[$i]=i;
   };
}
#Now store the tracts from the S' score file in a hash, using the alleles
# of the putatively archaic haplotype as the keys, and the tract ID as the
# value:
filenum==1&&FNR>1{
   tract[$sprimecols["CHROM"]":"$sprimecols["POS"],$sprimecols["ALLELE"]]=$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"];
   if (length(debug) > 0) {
      print $sprimecols["CHROM"]":"$sprimecols["POS"]","$sprimecols["ALLELE"]" => "$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"] > "/dev/stderr";
   };
}
#The second file is the output of bcftools query -H -f '%CHROM:%POS[\t%GT]\n'
# so we trim off the prefixed "# " of the header line, and then removing
# the prefixed column numbers in order to get the column names:
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      sub("# ", "", $i);
      gsub("[[0-9]+]", "", $i);
      gsub(":GT", "", $i);
      querycols[i]=$i;
   };
}
#Determine the frequency of the archaic allele by first calculating all
# allele frequencies at the site, and then checking which of the alleles
# matches the archaic allele from S':
filenum==2&&FNR>1{
   AN=0;
   arcallele=-1;
   tractid="";
   for (i=2; i<=NF; i++) {
      ploidy=split($i, gt, "[/|]");
      for (j=1; j<=ploidy; j++) {
         AF[gt[j]]++;
         AN++;
         if (($1,gt[j]) in tract) {
            if (arcallele >= 0 && gt[j] != arcallele) {
               print "Somehow two archaic alleles were found at "$1": "arcallele" and "gt[j] > "/dev/stderr";
            };
            arcallele=gt[j];
            if (tractid != "" && tractid != tract[$1,gt[j]]) {
               print "Somehow two different tracts occurred at the same location: "$1" => "tractid" and "tract[$1,gt[j]] > "/dev/stderr";
            };
            tractid=tract[$1,gt[j]];
         };
      };
   };
   split($1, chrompos, ":");
   if (arcallele >= 0) {
      print chrompos[1], chrompos[2], tractid, AF[arcallele], AN, AF[arcallele]/AN, pop;
   } else {
#      print chrompos[1], chrompos[2], "NA", 0, AN, 0.0, pop;
   };
   delete AF;
}
