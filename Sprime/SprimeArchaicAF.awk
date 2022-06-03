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
# spop:   (required) The name of the populations the Sprime tracts came from
# pop:    (required) The name of the populations the genotypes belong to
# header: (optional) A flag indicating whether or not to print the header
#                    line (default: don't print header)
# allele: (optional) A flag indicating whether or not to print a column of
#                    the archaic allele (default: don't print arc allele)
# all:    (optional) A flag indicating whether or not to output all sites
#                    including those with archaic AF of 0
#                    (default: only print sites with arc AF > 0)
#The "header" flag facilitates combining the results across chromosomes and
# populations.
#Added the spop variable to enable simple concatenation of outputs across
# S' runs using different sample sets so that tracts are uniquely labeled.
BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   if (length(pop) == 0) {
      print "Missing pop variable, please set it." > "/dev/stderr";
      print "pop is the name of the query population" > "/dev/stderr";
      print "i.e. the source of the genotypes" > "/dev/stderr";
      exit 2;
   };
   if (length(spop) == 0) {
      print "Missing spop variable, please set it." > "/dev/stderr";
      print "spop is the name of the Sprime target population" > "/dev/stderr";
      print "i.e. the non-outgroup population passed to Sprime" > "/dev/stderr";
      exit 3;
   };
   if (length(header) > 0) {
      if (length(allele) > 0) {
         print "CHROM", "POS", "TractID", "ArchaicAllele", "ArchaicAlleleCount", "TotalAlleleCount", "ArchaicAlleleFrequency", "Population";
      } else {
         print "CHROM", "POS", "TractID", "ArchaicAlleleCount", "TotalAlleleCount", "ArchaicAlleleFrequency", "Population";
      };
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
   tract[$sprimecols["CHROM"]":"$sprimecols["POS"],$sprimecols["ALLELE"]]=spop"_"$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"];
   if (length(debug) > 0) {
      print $sprimecols["CHROM"]":"$sprimecols["POS"]","$sprimecols["ALLELE"]" => "spop"_"$sprimecols["CHROM"]"_"$sprimecols["SEGMENT"] > "/dev/stderr";
   };
   split($sprimecols["ALT"], sprimealleles, ",");
   sprimealleles[0]=$sprimecols["REF"];
   alleles[$sprimecols["CHROM"]":"$sprimecols["POS"]]=sprimealleles[$sprimecols["ALLELE"]+0];
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
         AC[gt[j]]++;
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
      if (length(allele) > 0) {
         print chrompos[1], chrompos[2], tractid, alleles[$1], AC[arcallele], AN, AC[arcallele]/AN, pop;
      } else {
         print chrompos[1], chrompos[2], tractid, AC[arcallele], AN, AC[arcallele]/AN, pop;
      };
   } else {
      if (length(all) > 0) {
         if (length(allele) > 0) {
            print chrompos[1], chrompos[2], "NA", "NA", 0, AN, 0.0, pop;
         } else {
            print chrompos[1], chrompos[2], "NA", 0, AN, 0.0, pop;
         };
      };
   };
   delete AC;
}
