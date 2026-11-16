#!/bin/awk -f
#This script takes a HapMap-formatted recombination rate map file and a
# list of chromosomes and positions of variants from the same chromosome
# and outputs a mapfile in the format required by SparsePainter (i.e.
# header line and then one line per variant with two columns: physical
# position of the variant and genetic map position of the variant).
#Genetic map position is calculated by linear interpolation as needed.
#Options:
# precision: Output the genetic map position with this many digits
#            (default: 6)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(precision) == 0) {
      precision=6;
   } else if (precision < 0) {
      print "Invalid precision specified, must be >= 0 digits after the decimal. Quitting." > "/dev/stderr";
      exit 2;
   };
   #Construct the format string for the output given the desired precision:
   outputfmt=sprintf("%%i%s%%.%if\n", OFS, precision);
   #Recombination rate scale:
   #Typically recombination rates are specified in units of cM/Mbp, so we
   # need to divide the physical distance in bp by this scaling factor to
   # convert to Mbp.
   ratescale=1000000;
   #Initialize the genetic map arrays:
   #This avoids issues with interpolation when the first recombination rate
   # map marker comes after the first query variant.
   gmappos[0]=0;
   gmaprate[0]=0;
   gmapdist[0]=0;
   #And initialize an index for when we're interpolating at query variants:
   gmap_interp_idx=0;
   #Keep track of which file we're on:
   filenum=0;
}
FNR==1{
   filenum+=1;
}
#First file is the HapMap-formatted recombination rate map file:
filenum==1&&FNR==1{
   #Check the header for necessary columns:
   for (i=1; i<=NF; i++) {
      gmapcols[$i]=i;
   };
   if (NF < 3 || NF > 4) {
      print "HapMap-formatted recombination rate map files should have 3 or 4 columns, but this one has "NF". Are you sure "FILE" is the right file? Quitting." > "/dev/stderr";
      exit 3;
   };
   if (!("Chromosome" in gmapcols && "Position(bp)" in gmapcols && "Rate(cM/Mb)" in gmapcols) || (NF == 4 && !("Map(cM)" in gmapcols))) {
      print "Recombination rate map file "FILE" is missing one or more of the expected columns: Chromosome, Position(bp), Rate(cM/Mb), and possibly Map(cM). Please check that you input the correct file. Quitting." > "/dev/stderr";
      exit 3;
   };

}
filenum==1&&FNR>1{
   gmap_idx+=1;
   if (gmap_idx == 1) {
      gmapchr=$gmapcols["Chromosome"];
   } else if ($gmapcols["Chromosome"] != gmapchr) {
      print "Recombination rate map file "FILE" contains more than one chromosome. Please provide a single-chromosome recombination rate map file. Quitting." > "/dev/stderr";
      exit 4;
   };
   gmappos[gmap_idx]=$gmapcols["Position(bp)"];
   gmaprate[gmap_idx]=$gmapcols["Rate(cM/Mb)"];
   #Now calculate the current gmapdist based on delta gmappos and previous gmaprate:
   gmapdist[gmap_idx]=gmapdist[gmap_idx-1]+(gmappos[gmap_idx]-gmappos[gmap_idx-1])*gmaprate[gmap_idx-1]/ratescale;
   #We could cross-check this against the "Map(cM)" column, if provided.
}
#Second file is the list of variant chromosomes and positions
filenum==2{
   if (FNR==1) {
      #Output the header for SparsePainter's mapfile format:
      print "position(base)", "genetic_distance(centiMorgan)";
   };
   #Check that the chromosome of the variant matches the map:
   if ($1 != gmapchr) {
      print "Variant provided ("$1":"$2") is on a different chromosome from the recombination rate map provided ("gmapchr"). Quitting." > "/dev/stderr";
      exit 5;
   };
   #Move along the gmap arrays until the element on or after the variant in physical space:
   while (gmappos[gmap_interp_idx] < $2 && gmap_interp_idx < gmap_idx) {
      gmap_interp_idx+=1;
   };
   if (gmappos[gmap_interp_idx] < $2 && gmap_interp_idx == gmap_idx) {
      #Catch the case where the variant is past the end of the recombination rate map:
      gmap_interp=gmapdist[gmap_interp_idx]+($2-gmappos[gmap_interp_idx])*gmaprate[gmap_interp_idx]/ratescale;
      printf outputfmt, $2, gmap_interp;
   } else {
      if (gmappos[gmap_interp_idx] == $2) {
         #Don't bother interpolating if the variant is on a map marker:
         printf outputfmt, $2, gmapdist[gmap_interp_idx];
      } else {
         #Interpolate from the previous map marker as needed:
         gmap_interp=gmapdist[gmap_interp_idx-1]+($2-gmappos[gmap_interp_idx-1])*gmaprate[gmap_interp_idx-1]/ratescale;
         printf outputfmt, $2, gmap_interp;
      };
   };
}
