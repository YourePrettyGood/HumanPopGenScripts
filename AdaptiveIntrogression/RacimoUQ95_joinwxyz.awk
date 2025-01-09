#!/bin/awk -f
#This script takes uncompressed VCFs of (1) the high-depth archaic hominins and
# (2) modern humans with accompanying AC and AN tags (and AA tags for the
# first VCF) in order to compile a tabular dataset for calculation of Racimo
# et al. 2017's U_{A,B,C,D}(w,x,y,z) and Q95_{A,B,C,D}(x,y,z) statistics.
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(usepop)==0 || usepop < 1) {
      usepop=0;
   };
   #Possible populations to exclude from AC and AN sums:
   if (length(exclude) > 0) {
      n_excl=split(exclude, excludearr, ",");
      for (i=1; i<=n_excl; i++) {
         excluded_pops[excludearr[i]]=i;
      };
   };
   #Iterate over arrays in a consistent order:
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Output the header the first time:
   header=1;
   #Keep track of which file we're on:
   filenum=0;
}
#Keep track of which file we're on:
FNR==1{
   filenum+=1;
}
#First file is the VCF of the archaics with AA,AC,AF,AN tags filled,
# so we store the frequency of all possible alleles at a SNP, the AA, and
# the AAF of DEN and NEA:
filenum==1&&!/^#/{
   #Construct the set of alleles:
   if ($5 != ".") {
      nalt=split($5, alleles, ",");
   };
   alleles[0]=$4;
   nalleles=nalt+1;
   #Initialize the archaic frequency array elements for this site:
   arcAF[$1,$2,"DEN","A"]=0.0;
   arcAF[$1,$2,"DEN","C"]=0.0;
   arcAF[$1,$2,"DEN","G"]=0.0;
   arcAF[$1,$2,"DEN","T"]=0.0;
   arcAF[$1,$2,"NEA","A"]=0.0;
   arcAF[$1,$2,"NEA","C"]=0.0;
   arcAF[$1,$2,"NEA","G"]=0.0;
   arcAF[$1,$2,"NEA","T"]=0.0;
   #Parse the INFO string:
   ninfo=split($8, info, ";");
   #Make sure we set a sentinel in case AA not found or non-ACGT:
   aa=".";
   #Keep track of whether AF tags are available, as fully homozygous ref
   # sites apparently don't get AF set, as there are no ALT alleles:
   DENAFset=0;
   NEAAFset=0;
   #Find and parse the INFO/AA, INFO/AF_DEN, and INFO/AF_NEA tags:
   for (i in info) {
      split(info[i], tagparts, "=");
      if (tagparts[1] == "AA" && tagparts[2] ~ /[ACGTacgt]/) {
         aa=tagparts[2];
      } else if (tagparts[1] == "AF_DEN") {
         #Store the DEN ALT allele frequencies if found:
         sumdenaf=0.0;
         ndenaf=split(tagparts[2], denafs, ",");
         for (i=1; i<=ndenaf; i++) {
            sumdenaf+=denafs[i];
            arcAF[$1,$2,"DEN",toupper(alleles[i])]=denafs[i];
         };
         #Calculate REF allele frequency for DEN:
         arcAF[$1,$2,"DEN",toupper(alleles[0])]=1.0-sumdenaf;
         #Mark that DENAF has been set:
         DENAFset=1;
      } else if (tagparts[1] == "AF_NEA") {
         #Store the NEA ALT allele frequencies if found:
         sumneaaf=0.0;
         nneaaf=split(tagparts[2], neaafs, ",");
         for (i=1; i<=nneaaf; i++) {
            sumneaaf+=neaafs[i];
            arcAF[$1,$2,"NEA",toupper(alleles[i])]=neaafs[i];
         };
         #Calculate REF allele frequency for NEA:
         arcAF[$1,$2,"NEA",toupper(alleles[0])]=1.0-sumneaaf;
         #Mark that NEAAF has been set:
         NEAAFset=1;
      };
   };
   #Fill in the REF allele frequency if no AF tag was provided:
   if (!DENAFset) {
      arcAF[$1,$2,"DEN",toupper($4)]=1.0;
   };
   if (!NEAAFset) {
      arcAF[$1,$2,"NEA",toupper($4)]=1.0;
   };
   #Store ancestral allele information if usable:
   if (aa ~ /[ACGTacgt]/) {
      AA[$1,$2]=aa;
      #Store the major allele for DEN and NEA:
      MAF_DEN=0.0;
      MAF_NEA=0.0;
      for (i=0; i<=nalleles; i++) {
         if (arcAF[$1,$2,"DEN",toupper(alleles[i])] > MAF_DEN) {
            MAF_DEN=arcAF[$1,$2,"DEN",toupper(alleles[i])];
            Major[$1,$2,"DEN"]=toupper(alleles[i]);
         };
         if (arcAF[$1,$2,"NEA",toupper(alleles[i])] > MAF_NEA) {
            MAF_NEA=arcAF[$1,$2,"NEA",toupper(alleles[i])];
            Major[$1,$2,"NEA"]=toupper(alleles[i]);
         };
      };
   };
}
#Second file is the VCF of modern humans with AC,AF,AN tags filled for each
# population:
filenum==2&&!/^#/{
   #Construct the set of alleles:
   if ($5 != ".") {
      nalt=split($5, alleles, ",");
   };
   alleles[0]=$4;
   nalleles=nalt+1;
   #Only process if AA is present for this site:
   if (($1,$2) in AA) {
      delete AC;
      delete AN;
      delete pops;
      #Parse the INFO string:
      ninfo=split($8, info, ";");
      for (i in info) {
         split(info[i], tagparts, "=");
         if (tagparts[1] ~ /^AC_|AN_/) {
            ntagnameparts=split(tagparts[1], tagnameparts, "_");
            #If we want to keep AFs by pop rather than by superpop,
            # toggle usepop
            if (usepop && ntagnameparts == 3) {
               pop=tagparts[1];
               sub(/^A[CN]_/, "", pop);
            } else {
               pop=tagnameparts[2];
            };
            #Exclude any pops if requested:
            longpopid=tagparts[1];
            sub(/^A[CN]_/, "", longpopid);
            if (!(longpopid in excluded_pops)) {
               pops[pop]+=1;
               #Split out the ACs by ALT allele:
               if (tagnameparts[1] == "AC") {
                  naltcounts=split(tagparts[2], acs, ",");
                  for (i=1; i<=naltcounts; i++) {
                     AC[pop,"S"]+=acs[i];
                     AC[pop,alleles[i]]+=acs[i];
                  };
               } else if (tagnameparts[1] == "AN") {
                  AN[pop]+=tagparts[2];
               };
            };
         };
      };
      #Output the header once:
      if (header > 0) {
         printf "CHROM\tPOS\tAA";
         printf "\tDEN_MAJ\tDEN_A\tDEN_C\tDEN_G\tDEN_T";
         printf "\tNEA_MAJ\tNEA_A\tNEA_C\tNEA_G\tNEA_T";
         for (p in pops) {
            printf "\t%s_AN\t%s_A\t%s_C\t%s_G\t%s_T", p, p, p, p, p;
         };
         printf "\n";
         header=0;
      };
      #Output the archaic portions of the output line:
      printf "%s\t%i\t%s", $1, $2, AA[$1,$2];
      printf "\t%s\t%f\t%f\t%f\t%f", Major[$1,$2,"DEN"], arcAF[$1,$2,"DEN","A"], arcAF[$1,$2,"DEN","C"], arcAF[$1,$2,"DEN","G"], arcAF[$1,$2,"DEN","T"];
      printf "\t%s\t%f\t%f\t%f\t%f", Major[$1,$2,"NEA"], arcAF[$1,$2,"NEA","A"], arcAF[$1,$2,"NEA","C"], arcAF[$1,$2,"NEA","G"], arcAF[$1,$2,"NEA","T"];
      #Fill out the REF allele counts for each pop and output the line:
      for (p in pops) {
         if (AN[p] > 0) {
            AC[p,$4]=AN[p]-AC[p,"S"];
            printf "\t%i\t%f\t%f\t%f\t%f", AN[p], AC[p,"A"]/AN[p], AC[p,"C"]/AN[p], AC[p,"G"]/AN[p], AC[p,"T"]/AN[p];
         } else {
            printf "\t0\t0\t0\t0\t0";
         };
      };
      printf "\n";
   };
}
