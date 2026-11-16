#!/bin/awk -f
#This script takes a lifted annotation GFF3 produced by LiftOn and
# summarizes characteristics of each lifted gene and transcript,
# producing one line per transcript per gene copy in the annotation.
#Output includes the total length of exons and CDSes (if applicable)
# of the lifted transcript, nucleotide and amino acid identity
# (if applicable) of the lifted transcript and its coding sequence,
# gene type, and gene name (which may or may not differ from the gene
# ID).
#Tested with LiftOn GFF3s using GENCODE r50 GRCh38 GFF3 as source
# annotation, which includes protein-coding genes, ncRNAs, a wide
# variety of other non-protein-coding gene types, and pseudogenes.
#Output includes a header and is sorted in order of gene name, gene ID,
# transcript ID, and then chromosome and position. This ensures that
# transcripts from the same gene copy are grouped together.
#Options:
# progress:    Print number of input lines processed every X lines
#              to STDERR
#              (default: 1000000)
# no_progress: Do not print progress
#              (default: 0, i.e. DO print progress)
BEGIN{
   FS="\t";
   OFS=FS;
   tx_count=0;
   sort_cmd="sort -k4,4V -k5,5V -k8,8V -k1,1V -k2,2n -k3,3n";
   if (length(progress) == 0) {
      progress=1000000;
   };
   if (length(no_progress) == 0) {
      no_progress=0;
   };
}
!/^#/{
   n=split($9, tags, ";");
   id="";
   parent="";
   gene_name="";
   gene_type="";
   prot_nuc_ident="";
   prot_aa_ident="";
   extra_copies="";
   for (i in tags) {
      split(tags[i], tagparts, "=");
      if (tagparts[1] == "ID") {
         id=tagparts[2];
      } else if (tagparts[1] == "extra_copy_number") {
         extra_copies=tagparts[2];
      } else if (tagparts[1] == "Parent") {
         parent=tagparts[2];
      } else if (tagparts[1] == "gene_name") {
         gene_name=tagparts[2];
      } else if (tagparts[1] == "gene_type") {
         gene_type=tagparts[2];
      } else if (tagparts[1] == "dna_identity") {
         prot_nuc_ident=tagparts[2];
      } else if (tagparts[1] == "protein_identity") {
         prot_aa_ident=tagparts[2];
      };
   };
   if ($3 == "gene") {
      gene_chrom[id]=$1;
      gene_start[id]=$4;
      gene_end[id]=$5;
      gene_names[id]=gene_name;
      gene_types[id]=gene_type;
      if (gene_name in gene_copies) {
         if (extra_copies > gene_copies[gene_name]) {
            gene_copies[gene_name]=extra_copies;
         };
      } else {
         gene_copies[gene_name]=0;
      };
   } else if ($3 == "transcript" || $3 == "mRNA") {
      tx_order[tx_count]=id;
      tx_count+=1;
      tx_gene_map[id]=parent;
      prot_aa_idents[id]=prot_aa_ident;
      prot_nuc_idents[id]=prot_nuc_ident;
   } else if ($3 == "exon") {
      exon_lens[parent]+=$5-$4+1;
   } else if ($3 == "CDS") {
      cds_lens[parent]+=$5-$4+1;
   };
   if (no_progress < 1 && NR%progress==0) {
      print "Processed "NR" lines" > "/dev/stderr";
   };
}
END{
   PROCINFO["sorted_in"]="@ind_num_asc";
   print "Chrom", "Start", "End", "GeneName", "GeneID", "GeneType", "NumExtraCopies", "TxID", "TotalExonLen", "TotalCDSLen", "ProtNucIdent", "ProtAAIdent";
   for (i=0; i<tx_count; i++) {
      tx_id=tx_order[i];
      gene_id=tx_gene_map[tx_id];
      if (tx_id in cds_lens) {
         cds_len=cds_lens[tx_id];
      } else {
         cds_len=0;
      };
      print gene_chrom[gene_id], gene_start[gene_id], gene_end[gene_id], gene_names[gene_id], gene_id, gene_types[gene_id], gene_copies[gene_names[gene_id]], tx_id, exon_lens[tx_id], cds_len, prot_nuc_idents[tx_id], prot_aa_idents[tx_id] | sort_cmd;
   };
}
