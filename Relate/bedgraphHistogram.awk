#!/bin/awk -f
#This script generates an empirical CDF of a BEDGRAPH file.
#The empirical CDF is calculated from the right (i.e.
# largest values have smallest CDF values).
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   #Keep track of the total interval length per BEDGRAPH value:
   c[$4]+=$3-$2;
   #and the total length of intervals:
   t+=$3-$2;
}
END{
   PROCINFO["sorted_in"]="@ind_num_desc";
   #Calculate the cumulative site counts with decreasing threshold:
   cumul=0;
   for (i in c) {
      cumul+=c[i];
      print i, c[i], c[i]/t, cumul/t;
   };
}
