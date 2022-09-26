#!/bin/awk -f
#This script combines the counts from each per-chromosome run of
# bcftools +trio-stats for each combination of filter and trio.
#The output is thus a genome-wide estimate of the counts and
# rates of different kinds of trio outcomes, such as Mendelian
# errors, homozygous Mendelian errors, recurrent Mendelian
# errors, and novel singletons.
#Two rate estimates are provided for each, normalized either by
# the total number of passing sites for that trio, or the
# number of passing sites that contain at least one ALT allele
# amongst the trio members.
#The latter denominator is of course smaller, but excludes
# sites where all trio members are homozygous REF, which should
# be non-error sites.
#Homozygous Mendelian errors are likely to be genotyping errors,
# as they consist of a child homozygous for an allele not present
# in the parents.
#Recurrent Mendelian errors are cases where the putative error
# allele is found in other samples outside the trio, so may not
# be a genotyping error, or may be an error in the parents.
#Novel singletons are alleles only present in the child, and
# thus are counted as Mendelian errors.
BEGIN{
   FS="\t";
   OFS=FS;
}
/^FLT/{
   all[$1,$2,$3,$4]+=$5;
   havealt[$1,$2,$3,$4]+=$6;
   me[$1,$2,$3,$4]+=$7;
   novelsing[$1,$2,$3,$4]+=$8;
   untrans[$1,$2,$3,$4]+=$9;
   trans[$1,$2,$3,$4]+=$10;
   ts[$1,$2,$3,$4]+=$11;
   tv[$1,$2,$3,$4]+=$12;
   homme[$1,$2,$3,$4]+=$14;
   recurme[$1,$2,$3,$4]+=$15;
}
END{
   #Output in lexicographical order of filter and trio:
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Output a header:
   print "Filter", "Child", "Father", "Mother", "MendelianErrors", "HomozygMendelianErrors", "NovelSingletons", "RecurrentMendelianErrors", "ValidSitesWithAlt", "ValidSites", "MErate1", "MErate2", "HomMErate1", "HomMErate2", "NovelSrate1", "NovelSrate2", "RecurMErate1", "RecurMErate2", "Ts", "Tv", "Ts/Tv", "UntransmittedSingletons", "TransmittedSingletons";
   #Output each filter-trio combination:
   for (i in all) {
      split(i, a, SUBSEP);
      print a[1], a[2], a[3], a[4], me[i], homme[i], novelsing[i], recurme[i], havealt[i], all[i], me[i]/havealt[i], me[i]/all[i], homme[i]/havealt[i], homme[i]/all[i], novelsing[i]/havealt[i], novelsing[i]/all[i], recurme[i]/havealt[i], recurme[i]/all[i], ts[i], tv[i], ts[i]/tv[i], untrans[i], trans[i];
   };
}
