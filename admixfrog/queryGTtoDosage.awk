#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   #For bcftools query format intended, the GTs start at column 5:
   gtstart=5;
}
{
   if (length(header) > 0 && NR==1) {
      #Input has a header, so we need to reformat it:
      #First thing to do is eliminate the prefixed "# "
      $0=substr($0, 3);
      #Now we eliminate the "[#]" prefixed to each column name:
      gsub(/\[[0-9]\]/, "");
      #And finally this is one rough way to get rid of the ":GT" for all
      # the GT columns:
      gsub(/:GT/, "");
      #Now we print out the header line as adjusted:
      print $0;
   } else {
      for (i=gtstart; i<=NF; i++) {
         ploidy=split($i, gt, "[/|]");
         dose=gt[1];
         if (dose != ".") {
            for (j=2; j<=ploidy; j++) {
               dose+=gt[j];
            };
         };
         $i=dose;
      };
      print $0;
   };
}
