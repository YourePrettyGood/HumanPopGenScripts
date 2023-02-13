#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(source) == 0 || source == "tomahawk") {
      chroma="CHROM_A";
      posa="POS_A";
      posb="POS_B";
      rsquared="R2";
      chromb="CHROM_B";
   } else {
      chroma="CHR";
      posa="POS1";
      posb="POS2";
      rsquared="R^2";
      chromb="CHR";
   };
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1{
   keep[$1,$2]=1;
}
filenum==2&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   print "CHROM_A", "POS_A", "POS_B", "Rsquared", "CHROM_B";
}
filenum==2&&FNR>1{
   if ((($cols[chroma],$cols[posa]) in keep) && (($cols[chromb],$cols[posb]) in keep)) {
      print $cols[chroma], $cols[posa], $cols[posb], $cols[rsquared], $cols[chromb];
   };
}
