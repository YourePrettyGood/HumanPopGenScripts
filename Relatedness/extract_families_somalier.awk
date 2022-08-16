#!/bin/awk -f
#This script takes the *.samples.tsv and *.pairs.tsv files output by
# somalier relate --infer and extracts trios and other interesting
# relatedness patterns for further pedigree inference. somalier's
# built-in pedigree inference seems to get confused in certain cases
# so we build around some of these edge cases.
#For instance, we re-call the sex of individuals based on X and Y
# depth relative to the mean depth using the following criteria:
# 1) Male if Y scaled depth > 0.3 && 0.3 < X scaled depth < 0.7
# 2) Female if Y scaled depth < 0.1 && X scaled depth > 0.7
# 3) Unknown otherwise
#Furthermore, we allow flexible thresholds such as:
# po_thresh:     The supremum value of IBS0/IBS2 allowed for
#                parent-offspring pairs
#                (default: 0.005)
# unrel_thresh:  The supremum value of relatedness allowed for
#                unrelated pairs
#                (default: 0.06)
# fd_inf_rel:    The infimum value of relatedness allowed for
#                first-degree pairs
#                (default: 0.38)
# fd_sup_rel:    The supremum value of relatedness allowed for
#                first-degree pairs
#                (default: 0.62)
# fs_inf_ibs:    The infimum value of IBS0/IBS2 allowed for
#                full-sib pairs
#                (default: 0.015)
# fs_sup_ibs:    The supremum value of IBS0/IBS2 allowed for
#                full-sib pairs
#                (default: 0.052)
BEGIN{
   FS="\t";
   OFS=FS;
   #Set default thresholds:
   if (length(po_thresh) == 0) {
      po_thresh=0.005;
   };
   if (length(unrel_thresh) == 0) {
      unrel_thresh=0.06;
   };
   if (length(fd_inf_rel) == 0) {
      fd_inf_rel=0.38;
   };
   if (length(fd_sup_rel) == 0) {
      fd_sup_rel=0.62;
   };
   if (length(fs_inf_ibs) == 0) {
      fs_inf_ibs=0.015;
   };
   if (length(fs_sup_ibs) == 0) {
      fs_sup_ibs=0.052;
   };
   #Keep track of the input file index:
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1&&FNR==1{
   sub("#", "", $0);
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum==1&&FNR>1{
   if ($cols["Y_depth_mean"]/$cols["depth_mean"] > 0.3 && $cols["X_depth_mean"]/$cols["depth_mean"] < 0.7 && $cols["X_depth_mean"]/$cols["depth_mean"] > 0.3) {
      sex[$cols["sample_id"]]="M";
   } else if ($cols["Y_depth_mean"]/$cols["depth_mean"] < 0.1 && $cols["X_depth_mean"]/$cols["depth_mean"] > 0.7) {
      sex[$cols["sample_id"]]="F";
   } else {
      sex[$cols["sample_id"]]="U";
   };
}
filenum==2&&FNR==1{
   delete cols;
   sub("#", "", $0);
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
filenum==2&&FNR>1{
   hommatch=$cols["ibs0"]/$cols["ibs2"];
   if ($cols["relatedness"] > fd_inf_rel && $cols["relatedness"] < fd_sup_rel) {
      if (hommatch < po_thresh) {
         if ($cols["sample_a"] < $cols["sample_b"]) {
            po[$cols["sample_a"],$cols["sample_b"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
         } else {
            po[$cols["sample_b"],$cols["sample_a"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
         };
      } else if (hommatch > fs_inf_ibs && hommatch < fs_sup_ibs) {
         fs[$cols["sample_a"],$cols["sample_b"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
         fs[$cols["sample_b"],$cols["sample_a"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
      };
   } else if ($cols["relatedness"] < unrel_thresh) {
      unrel[$cols["sample_a"],$cols["sample_b"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
      unrel[$cols["sample_b"],$cols["sample_a"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
   };
   all[$cols["sample_a"],$cols["sample_b"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
   all[$cols["sample_b"],$cols["sample_a"]]=sprintf("%f,%f", $cols["relatedness"], hommatch);
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   for (i in po) {
      split(i, a, SUBSEP);
      for (j in po) {
         split(j, b, SUBSEP);
         diff="";
         common="";
         if (i != j) {
            if (a[1] == b[1]) {
               diff=a[2] SUBSEP b[2];
               common=a[1];
            } else if (a[2] == b[1]) {
               diff=a[1] SUBSEP b[2];
               common=a[2];
            } else if (a[1] == b[2]) {
               diff=a[2] SUBSEP b[1];
               common=a[1];
            } else if (a[2] == b[2]) {
               diff=a[1] SUBSEP b[1];
               common=a[2];
            } else {
               continue;
            }
            if (diff in unrel) {
               split(diff, parents, SUBSEP);
               if (sex[parents[1]] == sex[parents[2]] || sex[parents[1]] == "U" || sex[parents[2]] == "U") {
                  print "Warning: Sexes of putative parents in trio "diff" don't fit expectations";
               };
               #Only store the maternal first version of each trio:
               if (sex[parents[1]] == "M") {
                  continue;
               };
               if (diff in trio) {
                  if (!((trio[diff],common) in fs)) {
                     print "Warning: Putative children of trio "diff" aren't full sibs";
                  };
                  trio[diff]=trio[diff]";"common;
               } else {
                  trio[diff]=common;
               };
            } else if (diff in fs) {
               split(diff, children, SUBSEP);
               #Only output one permutation of the set:
               if (children[1] < children[2]) {
                  print "SameParent", common, sex[common], children[1], children[2];
               };
            } else {
               split(diff, different, SUBSEP);
               #Only output one permutation of the set:
               if (different[1] < different[2]) {
                  print "Unknown", common, sex[common], different[1], sex[different[1]], different[2], sex[different[2]], all[diff], all[common,different[1]], all[common,different[2]];
               };
            };
         };
      };
   };
   for (i in trio) {
      split(i, a, SUBSEP);
      split(trio[i], b, ";");
      childstats="";
      for (j in b) {
         childstats=childstats OFS all[b[j],a[1]] OFS all[b[j],a[2]];
      };
      print "trio", a[1], sex[a[1]], a[2], sex[a[2]], trio[i], sex[trio[i]], all[i] childstats;
   };
   for (i in po) {
      split(i, a, SUBSEP);
      split(po[i], b, ",");
      print a[1], sex[a[1]], a[2], sex[a[2]], "PO", b[1], b[2];
   };
   for (j in fs) {
      split(j, a, SUBSEP);
      if (a[1] < a[2]) {
         split(fs[j], b, SUBSEP);
         print a[1], sex[a[1]], a[2], sex[a[2]], "FS", b[1], b[2];
      };
   };
}
