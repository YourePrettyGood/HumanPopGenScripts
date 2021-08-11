#!/bin/awk -f

BEGIN{
   FS="\t";
   OFS=FS;
   filenum=0;
   PROCINFO["sorted_in"]="@ind_num_asc";
}
FNR==1{
   filenum++;
}
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      metacols[$i]=i;
   };
}
filenum==1&&FNR>1{
   region[$metacols["Sample"]]=$metacols["Region"];
}
filenum==2&&/^#CHROM/{
   for (i=10; i<=NF; i++) {
      if (region[$i] == "Denisovan" || region[$i] == "Neandertal") {
         useindiv[i]=region[$i];
         indiv[i]=$i;
      };
   };
}
filenum==2&&!/^#/{
   n_alts=split($5, alleles, ",");
   alleles[0]=$4;
   sitekey=$1":"$2;
   split($9, format, ":");
   for (f in format) {
      if (format[f] == "GT") {
         gtindex=f;
         break;
      };
   };
   for (i in useindiv) {
      split($i, samplearr, ":");
      ploidy=split(samplearr[gtindex], gt, "[/|]");
      for (j=1; j<=ploidy; j++) {
         if (gt[j] == ".") {
            continue;
         };
         if (!((useindiv[i],alleles[gt[j]]) in AC)) {
            if ((sitekey,useindiv[i]) in arcalleles) {
               arcalleles[sitekey,useindiv[i]]=arcalleles[sitekey,useindiv[i]]","alleles[gt[j]];
            } else {
               arcalleles[sitekey,useindiv[i]]=alleles[gt[j]];
            };
         };
         AC[useindiv[i],alleles[gt[j]]]++;
         if (!((indiv[i],alleles[gt[j]]) in AC)) {
            if ((sitekey,indiv[i]) in arcalleles) {
               arcalleles[sitekey,indiv[i]]=arcalleles[sitekey,indiv[i]]","alleles[gt[j]];
            } else {
               arcalleles[sitekey,indiv[i]]=alleles[gt[j]];
            };
         };
         AC[indiv[i],alleles[gt[j]]]++;
      };
   };
   delete AC;
}
filenum==3&&FNR==1{
   for (i=1; i<=NF; i++) {
      sprimecols[$i]=i;
   };
   print $0, "NeandertalMatch", "DenisovanMatch", "NeandertalAlleles", "DenisovanAlleles", "AltaiMatch", "VindijaMatch", "ChagyrskayaMatch";
}
filenum==3&&FNR>1{
   sitekey=$sprimecols["CHROM"]":"$sprimecols["POS"];
   split($sprimecols["ALT"], alleles, ",");
   alleles[0]=$sprimecols["REF"];
   sprimeallele=alleles[$sprimecols["ALLELE"]];
   Denisovan="missing";
   Neandertal="missing";
   Altai="missing";
   Vindija="missing";
   Chagyrskaya="missing";
   Denisovanallelestr="Unk";
   Neandertalallelestr="Unk";
   if ((sitekey,"Denisovan") in arcalleles) {
      Denisovan="mismatch";
      Denisovanallelestr=arcalleles[sitekey,"Denisovan"];
      split(arcalleles[sitekey,"Denisovan"], denalleles, ",");
      for (a in denalleles) {
         if (denalleles[a] == sprimeallele) {
            Denisovan="match";
         };
      };
   };
   if ((sitekey,"Neandertal") in arcalleles) {
      Neandertal="mismatch";
      Neandertalallelestr=arcalleles[sitekey,"Neandertal"];
      split(arcalleles[sitekey,"Neandertal"], neandalleles, ",");
      for (a in neandalleles) {
         if (neandalleles[a] == sprimeallele) {
            Neandertal="match";
         };
      };
   };
   if ((sitekey,"AltaiNeandertal") in arcalleles) {
      Altai="mismatch";
      split(arcalleles[sitekey,"AltaiNeandertal"], altaialleles, ",");
      for (a in altaialleles) {
         if (altaialleles[a] == sprimeallele) {
            Altai="match";
         };
      };
   };
   if ((sitekey,"Vindija33.19") in arcalleles) {
      Vindija="mismatch";
      split(arcalleles[sitekey,"Vindija33.19"], vindijaalleles, ",");
      for (a in vindijaalleles) {
         if (vindijaalleles[a] == sprimeallele) {
            Vindija="match";
         };
      };
   };
   if ((sitekey,"Chagyrskaya-Phalanx") in arcalleles) {
      Chagyrskaya="mismatch";
      split(arcalleles[sitekey,"Chagyrskaya-Phalanx"], chagyrskayaalleles, ",");
      for (a in chagyrskayaalleles) {
         if (chagyrskayaalleles[a] == sprimeallele) {
            Chagyrskaya="match";
         };
      };
   };
   print $0, Neandertal, Denisovan, Neandertalallelestr, Denisovanallelestr, Altai, Vindija, Chagyrskaya;
}
