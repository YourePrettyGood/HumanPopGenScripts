#!/bin/awk -f
#This script takes the output of summarize_lifted_genes.awk and summarizes
# the multiple transcripts of each gene down to a single row per gene
# indicating how many of the predicted extra copies of the gene in the
# target assembly have sufficiently complete lifted annotations.
#"Sufficiently complete" here means that all transcripts from the
# putatively orthologous lifted gene copy are also represented in the
# extra copy with total exon length (and total CDS length, if applicable)
# covering at least len_thresh (default 0.95) of the putatively orthologous
# copy's value, and sequence nucleotide identity (and amino acid identity,
# if applicable) are at least ident_thresh (default 0.9) of the putatively
# orthologous copy's value.
#Output can be further filtered to only include genes with at least
# min_copies (default 0) extra copies annotated, and at least min_full_copies
# (default 0) sufficiently complete extra copies.
#The defaults output all genes, even if they do not have extra copies
# annotated.
#Options:
# len_thresh:      Minimum proportion of total exon (and CDS, if applicable)
#                  length required for an extra copy to be considered full
#                  (default: 0.95)
# ident_thresh:    Minimum proportion of nucleotide (and amino acid, if applicable)
#                  identity required for an extra copy to be considered full
#                  (default: 0.9)
# min_copies:      Minimum number of annotated extra gene copies for a gene
#                  to be output
#                  (default: 0, i.e. output all genes, even if no extra copies)
# min_full_copies: Minimum number of annotated extra gene copies evaluated as
#                  being full in order for a gene to be output
#                  (default: 0, i.e. output genes even if their extra copies are
#                   only partial)
# asm:             Assembly name to include as a separate column
#                  (default: "")
BEGIN{
   FS="\t";
   OFS=FS;
   gene_count=0;
   if (length(len_thresh) == 0) {
      len_thresh=0.95;
   };
   if (length(ident_thresh) == 0) {
      ident_thresh=0.9;
   };
   if (length(min_copies) == 0) {
      min_copies=0;
   };
   if (length(min_full_copies) == 0) {
      min_full_copies=0;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   n=split($cols["GeneID"], idparts, "_");
   gene_id=idparts[1];
   n_tx_parts=split($cols["TxID"], txidparts, "_");
   tx_id=txidparts[1];
   if (n == 1) {
      gene_order[gene_id]=gene_count;
      gene_count+=1;
      ortho_stats[gene_id,"GeneName"]=$cols["GeneName"];
      ortho_stats[gene_id,"GeneType"]=$cols["GeneType"];
      ortho_stats[gene_id,tx_id,"TotalExonLen"]=$cols["TotalExonLen"];
      ortho_stats[gene_id,tx_id,"TotalCDSLen"]=$cols["TotalCDSLen"];
      ortho_stats[gene_id,tx_id,"ProtNucIdent"]=$cols["ProtNucIdent"];
      ortho_stats[gene_id,tx_id,"ProtAAIdent"]=$cols["ProtAAIdent"];
      ortho_stats[gene_id,"NumTx"]+=1;
      gene_tx_map[gene_id,ortho_stats[gene_id,"NumTx"]]=tx_id;
      ortho_stats[gene_id,"NumExtraCopies"]=$cols["NumExtraCopies"];
   } else {
      if (gene_id SUBSEP "GeneType" in ortho_stats) {
         if (ortho_stats[gene_id,tx_id,"TotalExonLen"] > 0) {
            total_exon_len_ratio=$cols["TotalExonLen"]/ortho_stats[gene_id,tx_id,"TotalExonLen"];
         } else {
            #For now, set to 0, although this can be caused by the putative ortholog
            # missing a transcript:
            total_exon_len_ratio=0;
         };
         if (ortho_stats[gene_id,tx_id,"TotalCDSLen"] > 0) {
            total_cds_len_ratio=$cols["TotalCDSLen"]/ortho_stats[gene_id,tx_id,"TotalCDSLen"];
         } else {
            total_cds_len_ratio="NA";
         };
         if (ortho_stats[gene_id,tx_id,"ProtNucIdent"] > 0) {
            prot_nuc_ident_ratio=$cols["ProtNucIdent"]/ortho_stats[gene_id,tx_id,"ProtNucIdent"];
         } else {
            #This generally shouldn't occur, but might occur due to the putative ortholog
            # missing a transcript:
            prot_nuc_ident_ratio=0;
         };
         if (ortho_stats[gene_id,tx_id,"ProtAAIdent"] > 0) {
            prot_aa_ident_ratio=$cols["ProtAAIdent"]/ortho_stats[gene_id,tx_id,"ProtAAIdent"];
         } else {
            prot_aa_ident_ratio="NA";
         };
         if (ortho_stats[gene_id,"GeneType"] == "protein_coding") {
            if (total_exon_len_ratio > len_thresh && total_cds_len_ratio > len_thresh) {
               ortho_stats[gene_id,idparts[2],"FullLen",tx_id]+=1;
            };
            if (prot_nuc_ident_ratio > ident_thresh && prot_aa_ident_ratio > ident_thresh) {
               ortho_stats[gene_id,idparts[2],"GoodIdent",tx_id]+=1;
            };
         } else {
            if (total_exon_len_ratio > len_thresh) {
               ortho_stats[gene_id,idparts[2],"FullLen",tx_id]+=1;
            };
            if (prot_nuc_ident_ratio > ident_thresh) {
               ortho_stats[gene_id,idparts[2],"GoodIdent",tx_id]+=1;
            };
         };
      } else {
         print "Lifted gene copy appears to be listed before the putatively orthologous lifted gene for "idparts[1]" ("$cols["GeneName"]", "$cols["GeneID"]"), this is unexpected." > "/dev/stderr";
         exit 3;
      };
   };
}
END{
   PROCINFO["sorted_in"]="@val_num_asc";
   print "GeneName", "GeneID", "GeneType", "NumTx", "NumExtraCopies", "NumFullCopies", "AsmName";
   for (gene_id in gene_order) {
      n_full_copies=0;
      n_copies=ortho_stats[gene_id,"NumExtraCopies"];
      for (i=1; i<=n_copies; i++) {
         full_tx=0;
         n_tx=ortho_stats[gene_id,"NumTx"];
         for (j=1; j<=n_tx; j++) {
            tx_id=gene_tx_map[gene_id,j];
#            print gene_id, tx_id, ortho_stats[gene_id,i,"FullLen",tx_id], ortho_stats[gene_id,i,"GoodIdent",tx_id];
            if (ortho_stats[gene_id,i,"FullLen",tx_id] > 0 && ortho_stats[gene_id,i,"GoodIdent",tx_id] > 0) {
               full_tx+=1;
            };
         };
         if (full_tx == n_tx) {
            n_full_copies+=1;
         };
      };
      if (n_copies >= min_copies && n_full_copies >= min_full_copies) {
         print ortho_stats[gene_id,"GeneName"], gene_id, ortho_stats[gene_id,"GeneType"], ortho_stats[gene_id,"NumTx"], ortho_stats[gene_id,"NumExtraCopies"], n_full_copies, asm;
      };
   };
}
