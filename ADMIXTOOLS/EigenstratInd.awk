#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(idcol) == 0) {
      print "Missing idcol variable, please specify it." > "/dev/stderr";
      exit 2;
   };
   if (length(metacol) == 0) {
      print "Missing metacol variable, please specify it." > "/dev/stderr";
      exit 3;
   };
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1&&FNR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum==1&&FNR>1{
   region[$cols[idcol]]=$cols[metacol];
}
filenum==2{
   if ($1 in region) {
      print $1, "U", region[$1];
   } else {
      print $1, "U", "9";
   };
}
