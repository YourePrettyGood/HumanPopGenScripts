#!/bin/awk -f
#This script summarizes the output of combineCrossCoal.py by performing
# the conversions from scaled time intervals and coalescence rates to
# generations, years, N_e, and relative cross-coalescence rate (rCCR).
#The input for this script is the output of combineCrossCoal.py.
#Optional arguments:
# gen_time:  Generation time (in years)
#            (default=29)
# mu:        Mutation rate (germline, per generation per basepair)
#            (default=0.0000000125)
#            (default from Scally & Durbin 2012 NRG doi: 10.1038/nrg3295
# print_gen: Include columns with time intervals in generations
#            (default=Don't include)
# Ne_eqn:    Use the original Ne equations from the MSMC2 guide rather
#            than the rearranged ones
#            (default=Use rearranged equations)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(gen_time) == 0 || gen_time <= 0) {
      gen_time=29;
   };
   if (length(mu) == 0 || mu <= 0) {
      mu=0.0000000125;
   };
}
#Keep track of the column names and print a header:
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Print the header line:
   if (length(print_gen) > 0 && print_gen > 0) {
      print "t_left_gens", "t_right_gens", "t_left", "t_right", "qNe", "tNe", "rCCR";
   } else {
      print "t_left", "t_right", "qNe", "tNe", "rCCR";
   };
}
#Perform the unit conversions and print out each converted line:
NR>1{
   #Convert the scaled time units to generations:
   t_left=$cols["left_time_boundary"]/mu;
   t_right=$cols["right_time_boundary"]/mu;
   #Calculate Ne using the guide.md way or the rearranged way:
   if (length(Ne_eqn) > 0 && Ne_eqn > 0) {
      qNe=(1/$cols["lambda_00"])/(2*mu);
      tNe=(1/$cols["lambda_11"])/(2*mu);
   } else {
      qNe=1/($cols["lambda_00"]*2*mu);
      tNe=1/($cols["lambda_11"]*2*mu);
   };
   #Calculate the relative cross-coalescence rate (rCCR):
   rCCR=(2*$cols["lambda_01"])/($cols["lambda_00"]+$cols["lambda_11"]);
   #Print out the converted values:
   if (length(print_gen) > 0 && print_gen > 0) {
      print t_left, t_right, t_left*gen_time, t_right*gen_time, qNe, tNe, rCCR;
   } else {
      print t_left*gen_time, t_right*gen_time, qNe, tNe, rCCR;
   };
}
