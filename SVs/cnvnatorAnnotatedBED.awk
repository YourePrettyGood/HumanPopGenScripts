#!/bin/awk -f
#This script does a quick bit of reformatting of the STDOUT
# from cnvnator -call into BED6 format where the Name column
# is analogous to the INFO column of a VCF.
#This output format enables easy processing with BEDtools.
#Arguments:
# ID: Sample ID, which isn't in the STDOUT but may be in
#     the filename
#Output annotations:
# SVTYPE: Uses the caps version of the first three letters
#         of column 1 of the input, usually DEL or DUP
# SVLEN:  Length of the called window in bp
# DPNORM: Normalized read depth across the window, 1 means
#         the adjusted average (diploid?) depth
# MQ0:    The fraction of reads in the window with MQ==0
# SNAME:  The sample ID, so what you provide with ID
BEGIN{
   FS="\t";
   OFS=FS;
}
NF==9{
   svtype=toupper(substr($1,1,3));
   split($2, chrompos, ":");
   split(chrompos[2], pos, "-");
   svlen=$3;
   normRD=$4;
   fracMQ0=$9;
   print chrompos[1], pos[1]-1, pos[2], "SVTYPE="svtype";SVLEN="svlen";DPNORM="normRD";MQ0="fracMQ0";SNAME="ID, ".", "+";
}
