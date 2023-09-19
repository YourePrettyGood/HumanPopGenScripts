#!/bin/awk -f
#
BEGIN{
   FS="\t";
   OFS=" ";
   if (length(idcol) == 0) {
      idcol="SampleID";
   };
   if (length(traitcol) == 0) {
      traitcol="Region";
   };
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1&&/^#/{
   for (i=2; i<=NF; i++) {
      sub("[[0-9]+]", "", $i);
      sub(":GT", "", $i);
      id[i]=$i;
   };
}
filenum==1&&!/^#/{
   ntaxa=0;
   for (i=2; i<=NF; i++) {
      ploidy=split($i, tgt, "[/|]");
      for (j=1; j<=ploidy; j++) {
         if (tgt[j] == ".") {
            seq[id[i],j]=seq[id[i],j]"N";
         } else {
            seq[id[i],j]=seq[id[i],j] tgt[j];
         };
         ntaxa++;
      };
   };
}
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum==2&&FNR>1{
   trait[$cols[idcol]]=$cols[traitcol];
   traitlist[$cols[traitcol]]++;
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Output the NEXUS header:
   print "#NEXUS";
   #Output the TAXA block:
   print "BEGIN TAXA;";
   print "DIMENSIONS NTAX="ntaxa";";
   printf "\n";
   print "TAXLABELS";
   for (i in seq) {
      split(i, a, SUBSEP);
      print a[1]"_"a[2];
      seqlen=length(seq[i]);
   };
   print ";";
   printf "\n";
   print "END;";
   printf "\n";
   #Output the CHARACTERS block:
   print "BEGIN CHARACTERS;";
   print "DIMENSIONS NCHAR="seqlen";";
   print "FORMAT DATATYPE=DNA MISSING=N GAP=- ;";
   print "MATRIX";
   printf "\n";
   for (i in seq) {
      split(i, a, SUBSEP);
      print a[1]"_"a[2], seq[i];
   };
   print ";";
   printf "\n";
   print "END;";
   printf "\n";
   #Output the TRAITS block:
   print "BEGIN TRAITS;";
   ntraits=0;
   traitstr="";
   for (t in traitlist) {
      ntraits++;
      traitstr=traitstr" "t;
      traitindex[t]=ntraits;
   };
   for (t in traitlist) {
      for (i=1; i<=ntraits; i++) {
         if (i == traitindex[t]) {
            traitcountstr[t]=traitcountstr[t]" 1";
         } else {
            traitcountstr[t]=traitcountstr[t]" 0";
         };
      };
   };
   print "  Dimensions NTRAITS="ntraits";";
   print "  Format labels=yes missing=? separator=spaces;";
   print "  TraitLabels"traitstr";";
   print "  Matrix";
   printf "\n";
   for (i in seq) {
      split(i, a, SUBSEP);
      print a[1]"_"a[2]""traitcountstr[trait[a[1]]];
   };
   print ";";
   printf "\n";
   print "END;";
}
