#!/bin/awk -f
#This script calculates population-specific allele counts given a
# metadata file and a VCF.
#The idea is that these frequencies are a precursor to identifying
# shared and population-private variants.
#
#The first input file should be a metadata file containing at a minimum
# columns for sample ID and population ID.
#Specify the names of these columns using the `idcol` and `metacol`
# arguments, respectively.
#The second input file is the VCF you want to process.
#
#Output is similar to the first five columns of a VCF, then a sequence
# of AC and AN columns for each population. The AC columns are
# populated by a comma-separated list of counts of each allele
# in the order they appear (i.e. REF,ALT[1],ALT[2],ALT[3],...)
#
#The `missinghomref` option treats missing genotypes as homozygous
# reference calls -- not ideal, but it's a somewhat common assumption.
#
#Changelog:
#Added an option to coerce missing genotypes to homozygous reference
BEGIN{
   FS="\t";
   OFS=FS;
   PROCINFO["sorted_in"]="@ind_str_asc";
   if (length(idcol) == 0) {
      print "Missing idcol variable, please set this to the column name for sample IDs in the first (i.e. metadata) file." > "/dev/stderr";
      exit 1;
   };
   if (length(metacol) == 0) {
      print "Missing metacol variable, please set this to the column name for the metadata by which you want samples partitioned." > "/dev/stderr";
      exit 2;
   };
   if (length(missinghomref) > 0) {
      print "Coercing missing genotypes to homozygous reference for counts." > "/dev/stderr";
   };
}
#First file is the metadata TSV with header
#Determine the column indices for sample IDs and the desired metadata:
FNR==NR&&FNR==1{
   for (i=1; i<=NF; i++) {
      if ($i == idcol) {
         id=i;
      } else if ($i == metacol) {
         meta=i;
      };
   };
   if (length(id) == 0) {
      print "No metadata column names matched "idcol", cannot proceed." > "/dev/stderr";
      exit 3;
   };
   if (length(meta) == 0) {
      print "No metadata column names matched "metacol", cannot proceed." > "/dev/stderr";
      exit 4;
   };
}
#Now we set up the samples hash to get sample IDs from column indices:
FNR==NR&&FNR>1{
   metadata[$id]=$meta;
}
#Proceed to the second file, which is a VCF
#We process the #CHROM header line to get sample order and
# check in case there are samples in the VCF missing from
# the metadata file:
FNR<NR&&/^#CHROM/{
   for (i=10; i<=NF; i++) {
      samples[i]=$i;
      if (length(metadata[$i]) == 0) {
         print "Sample "$i" not found in metadata file, did you input the correct VCF?" > "/dev/stderr";
         exit 5;
      };
      metadatavals[metadata[$i]]=1;
   };
   printf "%s\t%s\t%s\t%s\t%s", substr($1, 2), $2, $3, $4, $5;
   for (m in metadatavals) {
      printf "\tAC%s\tAN%s", m, m;
   };
   printf "\n";
}
#Now we process each VCF record and output the allele counts:
FNR<NR&&!/^#/{
#   AN=0;
   #Collect the alleles:
   n_alts=split($5, alleles, ",");
   alleles[0]=$4;
   #Optional detail for later -- parsing ancestral allele for polarization:
   if (length(aatag) > 0) {
      split($8, info, ";");
      for (i in info) {
         split(info[i], tag, "=");
         if (tag[1] == aatag) {
            AA=tag[2];
         };
      };
      #Deal with polarization here:
   };
   #Now figure out which FORMAT field is for GT:
   split($9, format, ":");
   for (i in format) {
      if (format[i] == "GT") {
         gtindex=i;
      };
   };
   #Iterate through the samples to fill out the allele counts:
   for (i=10; i<=NF; i++) {
      split($i, fields, ":");
      ploidy=split(fields[gtindex], gt, "[|/]");
      for (j=1; j<=ploidy; j++) {
         #In the case we want to coerce missing genotypes to homozygous ref:
         if (length(missinghomref) > 0 && gt[j] == ".") {
            gt[j]="0";
         };
         #This is a quick little trick to skip missing alleles:
         if (gt[j] in alleles) {
#            #Accumulate the global AC and AN counts:
#            AC[alleles[gt[j]]]++;
#            AN++;
            #Now accumulate the population-specific AC and AN:
            ACpop[metadata[samples[i]],alleles[gt[j]]]++;
            ANpop[metadata[samples[i]]]++;
         };
      };
   };
   #If we're simply outputting the population-specific AC and ANs,
   # use this section:
   #Print the feedthrough parts, so CHROM, POS, ID, REF, ALT:
   printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5;
   #Now the population-specific AC and ANs:
   for (p in metadatavals) {
      ACstr="";
      for (i=1; i<=n_alts; i++) {
         if (length(ACstr) > 0) {
            ACstr=ACstr",";
         };
         if ((p,alleles[i]) in ACpop) {
            ACstr=ACstr ACpop[p,alleles[i]];
         } else {
            ACstr=ACstr"0";
         };
      };
      if (p in ANpop) {
         printf "\t%s\t%i", ACstr, ANpop[p];
      } else {
         printf "\t%s\t0", ACstr;
      };
   };
   printf "\n";
#   delete AC;
   delete ACpop;
   delete ANpop;
}
