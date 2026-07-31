#!/bin/awk -f
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(src) == 0) {
      src="ABC";
      print "No enhancer src ID provided, defaulting to ABC gene-enhancer link column indices." > "/dev/stderr";
   };
   if (src != "ABC" && src != "Roadmap" && src != "EpiMap") {
      print "Invalid enhancer src ID provided. Quitting." > "/dev/stderr";
      exit 2;
   };
   if (src == "ABC") {
      tgt=8;
      cell=10;
      loc=11;
   } else if (src == "Roadmap") {
      tgt=8;
      cell=11;
      loc=12;
   } else if (src == "EpiMap") {
      tgt=8;
      cell=10;
      loc="";
   };
   if (length(outfmt) == 0) {
      outfmt="BED";
      print "No output format (outfmt) specified, defaulting to BED." > "/dev/stderr";
   };
   if (toupper(outfmt) != "BED" && toupper(outfmt) != "TSV") {
      print "Invalid output format (outfmt), must be either BED or TSV. Quitting." > "/dev/stderr";
      exit 3;
   };
   previd="";
   PROCINFO["sorted_in"]="@ind_str_asc";
}
{
   id=$1 OFS $2 OFS $3 OFS $4;
   if (id != previd) {
      if (previd != "") {
         if (toupper(outfmt) == "TSV") {
            n=split(previd, a, OFS);
            m=split(a[4], b, ";");
            svtype=b[1];
            sub("SVTYPE=", "", svtype);
            svlen=b[2];
            sub("SVLEN=", "", svlen);
            prefix=a[1] OFS a[2]+1 OFS svtype OFS svlen;
            targets="";
            for (t in targetarr) {
               targets=targets","t;
            };
            celltypes="";
            for (c in celltypearr) {
               celltypes=celltypes","c;
            };
            locs="";
            for (l in locarr) {
               locs=locs","l;
            };
         } else {
            prefix=previd;
         };
         suffix=substr(targets, 2) OFS substr(celltypes, 2);
         if (loc != "") {
            suffix=suffix OFS substr(locs, 2);
         };
         print prefix, suffix;
      };
      targets="";
      delete targetarr;
      celltypes="";
      delete celltypearr;
      locs="";
      delete locarr;
   };
   targets=targets","$tgt;
   targetarr[$tgt]+=1;
   celltypes=celltypes","$cell;
   celltypearr[$cell]+=1;
   if (loc != "") {
      locs=locs","$loc;
      locarr[$loc]+=1;
   };
   previd=id;
}
END{
   if (previd != "") {
      if (toupper(outfmt) == "TSV") {
         n=split(previd, a, OFS);
         m=split(a[4], b, ";");
         svtype=b[1];
         sub("SVTYPE=", "", svtype);
         svlen=b[2];
         sub("SVLEN=", "", svlen);
         prefix=a[1] OFS a[2]+1 OFS svtype OFS svlen;
         targets="";
         for (t in targetarr) {
            targets=targets","t;
         };
         celltypes="";
         for (c in celltypearr) {
            celltypes=celltypes","c;
         };
         locs="";
         for (l in locarr) {
            locs=locs","l;
         };
      } else {
         prefix=previd;
      };
      suffix=substr(targets, 2) OFS substr(celltypes, 2);
      if (loc != "") {
         suffix=suffix OFS substr(locs, 2);
      };
      print prefix, suffix;
   };
}
