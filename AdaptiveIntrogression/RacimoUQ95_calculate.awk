#!/bin/awk -f
#This script takes the output of RacimoUQ95_joinwxyz.awk and calculates
# the U_{A,B,C,D}(w,x,y,z) and Q95_{A,B,C,D}(w,y,z) statistics in windows
# based on the frequencies input. The key input parameters are A, B, C, D,
# w, x, y, z, and the window size. A and B can be comma-separated lists of
# populations (or superpopulations), while C and D default to NEA and DEN,
# respectively, though can be otherwise specified. A is the outgroup and B
# is the target, so w serves as a maximum/supremum on focal allele frequency
# in the outgroup, while x is the minimum/infimum of focal allele frequency
# in the target, although x is only used for the U statistic.
# y and z are exact frequency thresholds for C and D, respectively, so
# (y,z)=(1.0,0.0) corresponds to NEA-specific adaptive introgression by
# default, whereas (y,z)=(0.0,1.0) corresponds to DEN-specific.
#The output is a TSV of chromosome, window start and end (1-based),
# window number, the outgroup A, the target B, the number of assessed sites
# in the window, the U statistic, the Q95 statistic, and optionally a label
# indicating the archaic origin corresponding to the selection of (y,z).
BEGIN{
   FS="\t";
   OFS=FS;
   #Check and parse inputs:
   if (length(A) == 0) {
      print "A is missing, please specify a comma-separated list of pops you want to use as the outgroup A. Quitting." > "/dev/stderr";
      exit 2;
   };
   nA=split(A, Aarr, ",");
   for (i=1; i<=nA; i++) {
      Alist[Aarr[i]]=i;
   };
   if (length(B) == 0) {
      print "B is missing, please specify a single population to use as target B. Quitting." > "/dev/stderr";
      exit 3;
   };
   nB=split(B, Barr, ",");
   for (i=1; i<=nB; i++) {
      Blist[Barr[i]]=i;
      if (Barr[i] in Alist) {
         print "B ("Barr[i]") overlaps A ("A"), did you misspecify A? Quitting." > "/dev/stderr";
         exit 2;
      };
   };
   if (length(C) == 0) {
      print "C is missing, assuming you want NEA as C." > "/dev/stderr";
      C="NEA";
   };
   if (length(D) == 0) {
      print "D is missing, assuming you want DEN as D." > "/dev/stderr";
      D="DEN";
   };
   if (length(w) == 0) {
      print "w is missing, please specify the upper focal/derived allele frequency threshold in outgroup A ("A"). Quitting." > "/dev/stderr";
      exit 4;
   };
   if (w < 0 || w > 100) {
      print "w ("w") is outside of expected bounds, must be either in [0,1] or [0,100]. Quitting." > "/dev/stderr";
      exit 4;
   } else if (w > 1.0) {
      print "Interpreting w ("w") as a percentage and converting accordingly." > "/dev/stderr";
      w=w/100.0;
   };
   if (length(x) == 0) {
      print "x is missing, please specify the lower focal/derived allele frequency threshold in target B ("B"). Quitting." > "/dev/stderr";
      exit 5;
   };
   if (x < 0 || x > 100) {
      print "x ("x") is outside of expected bounds, must be either in [0,1] or [0,100]. Quitting." > "/dev/stderr";
      exit 5;
   } else if (x > 1.0) {
      print "Interpreting x ("x") as a percentage and converting accordingly." > "/dev/stderr";
      x=x/100.0;
   };
   if (length(y) == 0) {
      print "y is missing, please specify the focal/derived allele frequency in archaic C ("C"). Quitting." > "/dev/stderr";
      exit 6;
   };
   if (y < 0 || y > 100) {
      print "y ("y") is outside of expected bounds, must be either in [0,1] or [0,100]. Quitting" > "/dev/stderr";
      exit 6;
   } else if (y > 1.0) {
      print "Interpreting y ("y") as a percentage and converting accordingly." > "/dev/stderr";
      y=y/100.0;
   } else if (y != 0.0 && y != 1.0 && y != 100.0) {
      print "y is typically set to 0 or 1, yet you set it to "y". This is untested (mainly due to derived allele detection), so proceed at your own peril." > "/dev/stderr";
   };
   if (length(z) == 0) {
      print "z is missing, please specify the focal/derived allele frequency in archaic D ("D"). Quitting." > "/dev/stderr";
      exit 7;
   };
   if (z < 0 || z > 100) {
      print "z ("z") is outside of expected bounds, must be either in [0,1] or [0,100]. Quitting." > "/dev/stderr";
      exit 7;
   } else if (z > 1.0) {
      print "Interpreting z ("z") as a percentage and converting accordingly." > "/dev/stderr";
      z=z/100.0;
   } else if (z != 0.0 && z != 1.0 && z != 100.0) {
      print "z is typically set to 0 or 1, yet you set it to "z". This is untested (mainly due to derived allele detection), so proceed at your own peril." > "/dev/stderr";
   };
   if (length(window_size) == 0) {
      print "No window_size provided, defaulting to Racimo et al. 2017 choice of 40 kbp." > "/dev/stderr";
      window_size=40000;
   };
   if (length(header) == 0) {
      header=1;
   };
   #Initialize U_ABCD and Q_ABCD:
   U_ABCD=0;
   split("", Q_ABCD);
}
#A quick function to calculate the quantile of the Q array
# based on the R quantile(,type=1) approach (i.e. inverse empirical CDF):
function quantile(a, p) {
   if (length(a) == 0) {
      return "NA";
   };
   n=asort(a, q, "@val_num_asc");
   qi=int(n*p);
   if ((n*p)>qi) {
      qi+=1;
   };
   return q[qi];
}
NR==1{
   #Do some checking of the input arguments against the input file
   # and keep track of important column indices:
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
      if ($i ~ /_AN$/) {
         pop=$i;
         sub(/_AN$/, "", pop);
         pops[pop]+=1;
         weightcol[pop]=i;
         #We assume for now that the AF columns are in order after the AN column:
         freqcol[pop,"A"]=i+1;
         freqcol[pop,"C"]=i+2;
         freqcol[pop,"G"]=i+3;
         freqcol[pop,"T"]=i+4;
      };
   };
   Amissing=0;
   for (i=1; i<=nA; i++) {
      if (!(Aarr[i] in pops)) {
         Amissing+=1;
      };
   };
   if (Amissing > 0) {
      print "Missing "Amissing" of "nA" pops from outgroup A ("A") in the header of the input file. Quitting." > "/dev/stderr";
      exit 2;
   };
   Bmissing=0;
   for (i=1; i<=nB; i++) {
      if (!(Barr[i] in pops)) {
         Bmissing+=1;
      };
   };
   if (Bmissing > 0) {
      print "Missing "Bmissing" of "nB" pops from Target B ("B") in the header of the input file. Quitting." > "/dev/stderr";
      exit 3;
   };
   if (!(C"_MAJ" in cols)) {
      print "Archaic C ("C") not found amongst the pops in the header of the input file. Quitting." > "/dev/stderr";
      exit 8;
   };
   if (!(D"_MAJ" in cols)) {
      print "Archaic D ("D") not found amongst the pops in the header of the input file. Quitting." > "/dev/stderr";
      exit 9;
   };
   if (header > 0) {
      if (length(arclabel) > 0) {
         print "CHROM", "start", "end", "window", "A", "B", "num_sites", "U_ABCD", "Q95_ABCD", "Origin";
      } else {
         print "CHROM", "start", "end", "window", "A", "B", "num_sites", "U_ABCD", "Q95_ABCD";
      };
   };
}
NR>1{
   #Keep track of windows when possible, and output window summaries as new
   # windows arise:
   new_window=int(($cols["POS"]-1)/window_size)+1;
   if (new_window != window && length(window) != 0) {
      Q95_ABCD=quantile(Q_ABCD, 0.95);
      if (length(arclabel) > 0) {
         print $cols["CHROM"], (window-1)*window_size+1, window*window_size, window, A, B, sites, U_ABCD, Q95_ABCD, arclabel;
      } else {
         print $cols["CHROM"], (window-1)*window_size+1, window*window_size, window, A, B, sites, U_ABCD, Q95_ABCD;
      };
      U_ABCD=0;
      split("", Q_ABCD);
      sites=0;
   };
   window=new_window;
   #This site count allows us to keep track of windows with U=0:
   sites+=1;
   #Select the focal allele based on C and D:
   #We make a few assumptions here:
   # 1) y and/or z is large enough that the derived allele is the major
   #    allele in C or D
   # 2) There may be multiple derived alleles, but we care about one
   # 3) In the "both"/"Ambiguous" case, the major alleles of C and D must
   #    match
   #Given these assumptions, we partition (y,z) space based on the y=z line
   # and restricted to the upper-right quadrant (i.e. non-negative y and z),
   # and select the focal allele as the major allele of C if y >= z and the
   # major allele of D if y < z. We do not count the site if the focal allele
   # matches the ancestral allele or if both y and z are non-positive.
   #The idea is that in the generally useful case where y,z \in {0,1},
   # the major allele of the relevant archaic is the derived allele of
   # interest for the U and Q95 statistics.
   #If the archaics have the ancestral allele, the site isn't useful.
   #If the modern human haplotype of interest doesn't match the archaics,
   # there's something wonky going on, so skip the site, as it's probably
   # not adaptive introgression.
   focal_allele="";
   if (y > 0 && y >= z) {
      focal_allele=$cols[C"_MAJ"];
      if ($cols[C"_MAJ"] == $cols["AA"]) {
         next;
      };
      if (y == z && $cols[C"_MAJ"] != $cols[D"_MAJ"]) {
         next;
      };
   } else if (z > 0 && z > y) {
      focal_allele=$cols[D"_MAJ"];
      if ($cols[D"_MAJ"] == $cols["AA"]) {
         next;
      };
   } else {
      next;
   };
   #If the focal/derived allele is reasonable (i.e. A, C, G, or T),
   # calculate the frequencies in pops A-D and contribute to the U
   # count and the Q frequency array if frequency thresholds are
   # met:
   if (focal_allele ~ /^[ACGT]$/) {
      freqD=$cols[D"_"focal_allele];
      freqC=$cols[C"_"focal_allele];
      #freqB=$freqcol[B,focal_allele];
      #We do something a bit funky for A and B, since Racimo et al. 2017 allowed
      # for composite outgroups, so we have to calculate the frequency
      # of the focal allele in the composite outgroup as a weighted mean
      # across component populations where the weight is the AN of that
      # component population. This is algebraically equivalent to taking
      # the ratio of the sum of ACs to the sum of ANs, which is equal to
      # the AF of the composite outgroup if we handled it monolithically.
      AC_A=0;
      AN_A=0;
      for (i=1; i<=nA; i++) {
         AC_A+=$weightcol[Aarr[i]]*$freqcol[Aarr[i],focal_allele];
         AN_A+=$weightcol[Aarr[i]];
      };
      if (AN_A > 0) {
         freqA=AC_A/AN_A;
      } else {
         #Skip the site, since the frequency in A is unascertainable:
         next;
      };
      AC_B=0;
      AN_B=0;
      for (i=1; i<=nB; i++) {
         AC_B+=$weightcol[Barr[i]]*$freqcol[Barr[i],focal_allele];
         AN_B+=$weightcol[Barr[i]];
      };
      if (AN_B > 0) {
         freqB=AC_B/AN_B;
      } else {
         #Skip the site, since the frequency in B is unascertainable:
         next;
      };
      if (debug > 1 && freqA < w) {
         print "DEBUG:", $cols["CHROM"], $cols["POS"], $cols["AA"], focal_allele, freqA, freqB, freqC, freqD > "/dev/stderr";
      };
      #Filter primary on w, y, and z, since these apply to both U and Q95:
      if (freqD == z && freqC == y && freqA < w) {
         #Additionally filter on x for adding to the U count:
         if (freqB > x) {
            if (debug > 0 && debug < 2) {
               print $cols["CHROM"], $cols["POS"], $cols["AA"], focal_allele, freqA, freqB, freqC, freqD > "/dev/stderr";
            };
            U_ABCD+=1;
         };
         #Append freqB to the Q array:
         #Note: The index is arbitrary (though non-redundant), as quantile
         # will sort into a new sequentially-ordered array.
         Q_ABCD[$cols["POS"]]=freqB;
      };
   };
}
END{
   #Make sure we don't forget the last window:
   Q95_ABCD=quantile(Q_ABCD, 0.95);
   if (length(arclabel) > 0) {
      print $cols["CHROM"], (window-1)*window_size+1, window*window_size, window, A, B, sites, U_ABCD, Q95_ABCD, arclabel;
   } else {
      print $cols["CHROM"], (window-1)*window_size+1, window*window_size, window, A, B, sites, U_ABCD, Q95_ABCD;
   };
}
