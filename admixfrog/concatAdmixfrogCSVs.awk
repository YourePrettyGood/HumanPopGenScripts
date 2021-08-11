#!/bin/awk -f
BEGIN{
   FS=",";
   OFS=FS;
   filenum=0;
   PROCINFO["sorted_in"]="@ind_num_asc";
}
/^chrom/{
   filenum++;
}
filenum==1&&/^chrom/{
   for (i=1; i<=NF; i++) {
      basecols[i]=$i;
   };
   print;
}
filenum==1&&!/^chrom/{
   print;
}
filenum>1&&/^chrom/{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum>1&&!/^chrom/{
   for (i in basecols) {
      if (i>1) {
         printf "%s", OFS;
      };
      printf "%s", $cols[basecols[i]];
   };
   printf "\n";
}
