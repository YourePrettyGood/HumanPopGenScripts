#!/bin/awk -f
#This script processes the result of classification of a VCF to generate
# a rough BED of sites that are classified with a particular tag.
#One can use this to identify population-private variants, archaic
# population-private variants, singletons, etc.
#Required:
# targettag: Tag to select for
#Optional:
# incl_ref:  Also output sites where the tag applies for the REF
#            allele (generally undesirable, as the reference itself should
#            count when considering most tags)
BEGIN{
   FS="\t";
   OFS=FS;
   n_total_pops=0;
   if (length(targettag) == 0) {
      print "targettag not specified but is a required argument, quitting" > "/dev/stderr";
      exit 2;
   };
}
#Parse the header, looking for REF and tags, and also count the number of
# populations/regions in total based on the freqs_* columns:
NR==1{
   for (i=1; i<=NF; i++) {
      if ($i == "REF") {
         refcol=i;
      } else if ($i == "tags") {
         tagcol=i;
      } else if ($i ~ "freqs_") {
         n_total_pops++;
      };
   };
   if (length(refcol) == 0 || length(tagcol) == 0) {
      print "Couldn't find REF or tags column, quitting." > "/dev/stderr";
      exit 2;
   };
}
#Parse tags, split into (tagtype, pop(s), allele), ignore tags with REF allele,
# and count
NR>1{
   n_tags=split($tagcol, tags, ";");
   for (i=1; i<=n_tags; i++) {
      n_tag_parts=split(tags[i], tag_parts, "_");
      #Skip INVARIANT lines:
      if (n_tag_parts == 2 && tag_parts[1] == "INVARIANT") {
         next;
      };
      #Handle EXCLUDED tags here::
      if (tag_parts[1] == "EXCLUDED") {
         if (tags[i] == targettag) {
            print $1, $2-1, $2;
         };
         continue;
      };
      #All other tags should have at least 3 parts when separated by _:
      if (n_tag_parts < 3) {
         print "Site "$1":"$2" has an invalid tag with "n_tag_parts" parts: "tags[i] > "/dev/stderr";
         exit 3;
      };
      #Skip tags containing the reference allele so that we don't overinflate
      # the SHAREDALL class:
      if (tag_parts[n_tag_parts] == $refcol && length(incl_ref) == 0) {
         continue;
      };
      #Since the regions may contain _, we do loose tokenizing by
      # trimming the first and last element, and then split the poplist:
      poplist=tags[i];
      sub("^"tag_parts[1]"_", "", poplist);
      sub("_"tag_parts[n_tag_parts]"$", "", poplist);
      n_pops=split(poplist, pops, ",");
      tag_noallele=tags[i];
      sub("_"tag_parts[n_tag_parts]"$", "", tag_noallele);
      if (length(debug) > 0) {
         print tag_parts[1], poplist, tag_parts[n_tag_parts], tag_noallele, tags[i] > "/dev/stderr";
      };
      #If the current tag matches targettag, output a BED line:
#      if (n_pops == 1 && tag_noallele == targettag) {
      if (tag_noallele == targettag) {
         #Output CHROM\tPOS-1\tPOS\n for non-REF variants matching targettag:
         print $1, $2-1, $2;
      };
      if (length(debug) > 1) {
         print $1":"$2"_"i, tag_parts[1], poplist, tag_parts[n_tag_parts], n_tag_parts, n_pops, n_total_pops > "/dev/stderr";
      };
   };
}
