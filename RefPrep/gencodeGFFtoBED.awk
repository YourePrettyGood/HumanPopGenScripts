#!/bin/awk -f
#This script takes an uncompressed GFF3 annotation from GENCODE
# and outputs a BED file of intervals corresponding to gene
# records. The gene records may also be filtered by the gene_type
# tag, for instance to only retain protein-coding genes
# (gene_type == "protein_coding") or lncRNA genes
# (gene_type == "lncRNA").
#Optional arguments:
# types: Comma-separated list of gene_type values to retain
#        (default: protein_coding)
#        (if value is set to "all", then all gene_types are retained)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(types) == 0) {
      types="protein_coding";
   };
   if (types != "all") {
      split(types, typearr, ",");
      for (i in typearr) {
         gene_types[typearr[i]]=i;
      };
   };
   #sort command for the output:
   sort_cmd="sort -k1,1V -k2,2n -k3,3n";
}
#Only process gene records from the GFF3:
!/^#/&&$3=="gene"{
   #Split the GFF3 tag field into elements and extract the ID, gene_name, and gene_type:
   split($9, tags, ";");
   type="";
   id="";
   name="";
   for (i in tags) {
      split(tags[i], tagparts, "=");
      if (tagparts[1] == "ID") {
         #Get rid of the .# suffix of the ID for the gene:
         sub("[.][0-9]+$", "", tagparts[2]);
         id=tagparts[2];
      } else if (tagparts[1] == "gene_name") {
         name=tagparts[2];
      } else if (tagparts[1] == "gene_type") {
         type=tagparts[2];
      };
   };
   #Retain the selected gene_types:
   if (types == "all" || type in gene_types) {
      #Remember that GFF3 start is 1-based while BED start is 0-based, so subtract one.
      print $1, $4-1, $5, id, name | sort_cmd;
   };
}
