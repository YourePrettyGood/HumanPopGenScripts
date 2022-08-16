# HumanPopGenScripts
Scripts used for pre- and post-processing of various human population genomics analyses

## Types of analyses:

1. [admixfrog](#admixfrog)
1. [ADMIXTOOLS input prep](#admixtools)
1. [Archaic allele matching](#archaicmatch)
1. [De novo assembly](#assembly)
1. [Private allele counting](#privatealleles)
1. [Relatedness](#relatedness)
1. [Sprime](#sprime)
1. [SMC++](#smcpp)

<a name="admixfrog" />

## admixfrog

### `metadata_to_admixfrog_YAML.awk`

This script is fairly custom, but takes a population metadata file
and some optional arguments to define the parts of the admixfrog
YAML configuration file.

Arguments:

`idcol`: Column name for the sample ID column in the metadata file

`metacol`: Column name for the region column in the metadata file

`pseudohaploid`: (optional) Comma-separated list (no whitespace) of sample IDs that should be considered pseudohaploid

`neandertal`: (optional) Comma-separated list (no whitespace) of sample IDs to consider Neandertals, default below

`denisovan`: (optional) Comma-separated list (no whitespace) of sample IDs to consider Denisovans, default below

`outgroup`: (optional) Comma-separated list (no whitespace) of sample IDs to consider outgroups for polarization, default below

By default, Neandertals are automatically detected if their sample names
match "AltaiNeandertal", "Vindija33.19", or "Chagyrskaya-Phalanx", the
Denisovan is automatically detected if its sample name matches "Denisova", 
and the outgroup is "pan_troglodytes" by default (and is the only sample
assumed to be pseudohaploid). Recognized region codes that are mapped
to three-letter codes:

- Africa -> AFR
- America -> AMR
- CentralAsiaSiberia -> CAS
- Denisovan -> DEN
- EastAsia -> EAS
- Neandertal -> NEA
- Oceania_other -> OCE
- Oceania_PIB -> PIB
- Primates -> ANC
- SouthAsia -> SAS
- SoutheastAsia -> SEA
- WestEurasia -> EUR

### `renameBAMheader.awk`

This was just a silly script to make an hg19 BAM header compatible with
hs37d5 (as far as major chromosomes are concerned). I had to do this
because some samples were processed with these two slightly different
references.

Input file 1 is a 2-column TSV mapping hg19 scaffold names to hs37d5
scaffold names. Input file 2 is the SAM header for the BAM. You can
take the output of this script and pass it to `samtools reheader`.

### `concatAdmixfrogCSVs.awk`

This simple script lets you concatenate the CSVs produced by [admixfrog](https://github.com/BenjaminPeter/admixfrog).
The premise is that you can run admixfrog on each chromosome separately
to speed up/parallelize processing, and then combine the outputs as if
they were run all together. (Of course the inputs need to be decompressed,
and the output should be compressed.)

### `extractVariantSitesBED.awk`

This script does something super simple: Outputs BED3 lines of for
segregating sites in the input VCF.

At the moment, the last few lines assume diploid data, though the
code could easily be generalized to higher ploidy.

This script is mainly useful when you consider the archaic genomes
alone, and just look for sites segregating within the archaics.
Filtering down to only these sites reduces the noise introduced
to the admixfrog HMM. As per an e-mail with Ben Peter, this was
the procedure used for [Hajdinjak et al. 2021 Nature](https://doi.org/10.1038/s41586-021-03335-3).

### `simpleRLE.awk`

This script takes the posterior decoding from admixfrog (the .bin.xz file)
and performs a dirt-simple run-length encoding of the posterior decoded
state sequence. This can serve as an input to visualize the discretized
patterns of ancestry across chromosomes. The output includes one line per
tract of the same state and includes the coordinates of the tract in
three different coordinate spaces: SNP space ("id_*"), genetic map space
("map_*"), and physical/assembly space ("pos_*").

Usage:

```bash
xz -dc [path to admixfrog output].bin.xz | simpleRLE.awk | xz -c > [path to admixfrog output].crle.xz
```

<a name="admixtools" />

## ADMIXTOOLS

### `EigenstratInd.awk`

This script takes a tab-separated metadata file (with header) and a
sample ID list, and generates the .ind file for ADMIXTOOLS rather than
relying on `CONVERTF`.

Arguments:

`idcol`: Column name for the sample ID column in the metadata file

`metacol`: Column name for the region column in the metadata file

### `VCFtoEigenstrat.awk`

This script generates the .geno and .snp files for ADMIXTOOLS from
the output of `bcftools query -f '%CHROM:%POS\t%CHROM\t0.0\t%POS\t%REF\t%ALT\n[%GT\t]\n' [input VCF]`

Arguments:

`snp`: File name/path for the output .snp file

`geno`: File name/path for the output .geno file

`recratescale`: (optional) bp per cM to use for physical position scaling in the .snp file, default: 1000000 (i.e. 1 Mbp per cM)

<a name="archaicmatch" />

## Archaic allele matching


<a name="assembly" />

## De novo assembly

### `coordToBED.awk`

This script converts the output of [nucmer](https://github.com/mummer4/mummer)'s
`show-coords` program (with flags `-cdHloqT`) into a BED6 file in the
coordinate space of the reference used for `nucmer`. By default, the BED6
columns are:

1. Reference sequence name
1. 0-based alignment start in reference coordinates
1. 1-based alignment end in reference coordinates
1. GFF3-like tag string consisting of the query sequence name and overall length
1. Percent identity of the alignment
1. Orientation of the alignment

If you specify the `flip` argument (i.e. `-v "flip=1"`), then the coordinate
space is flipped to that of the query, hence the BED6 columns are:

1. Query sequence name
1. 0-based alignment start in query coordinates
1. 1-based alignment end in query coordinates
1. GFF3-like tag string consisting of the reference sequence name and overall length
1. Percent identity of the alignment
1. Orientation of the alignment

### `BEDbestHit.awk`

This script processes the BED produced by `coordToBED.awk` and sums alignment
lengths for each query-reference pair, then outputs a list of these pairs
sorted by their approximate total aligned length. The output includes an
estimate of the query coverage as well as the majority vote orientation.

Output columns:

1. Query sequence name
1. Reference sequence name
1. Approximate query coverage (can be > 100% due to a couple of reasons)
1. Sum of aligned lengths
1. Query sequence length
1. Majority vote orientation
1. Positive orientation votes
1. Negative orientation votes

For example, if you want to identify which query contigs belong to which
reference scaffolds, take a look at the output of:

```show-coords -cdHloqT [delta file] | coordToBED.awk | sort -k4,4V -k1,1V -k2,2n -k3,3n | BEDbestHit.awk```

You may want to re-sort the output by the query coverage estimate, but the
default sort order is nice because you can see if the longest query contigs
have any misassemblies. The query coverage estimate is probably also more
reliable when the aligned length is longer.

<a name="privatealleles" />

## Private allele counting

### `vcfToPerPopACAN.awk`

This script takes a population metadata file and a VCF (nominally of
biallelic variants, but may work with multiallelic) as input, and determines
the counts of alleles within each population (and calculates the
corresponding frequencies). Note that the counts and frequencies include
the reference allele first.

If you're post-processing the output, make sure the input doesn't include
invariant sites (e.g. `bcftools view -c 1 -m 2 [VCF]`). Invariant sites
don't cause errors in this script, but they're just a waste of time
to process most of the time.

Arguments:

`idcol`: Column name for the sample ID column in the metadata file

`metacol`: Column name for the region column in the metadata file

`missinghomref`: (optional) Coerce missing genotypes to count as homozygous reference in the allele counts and frequencies

### `bcftoolsQuerySiteClassification.awk`

This script takes the output of `vcfToPerPopACAN.awk` and classifies
each site based on sharing between populations, as well as optional
annotations for singletons and presence in a single outgroup population.
For instance, if you wanted to identify alleles not present in some
outgroup population "Outgroup", you could pass in `-v "exclude=Outgroup"`
and any alleles present in "Outgroup" would be marked with the tag
"EXCLUDED_Outgroup".

By default, the annotation tags include:

1. "PRIVATE_[population]_[allele]"
1. "SHARED_[population 1],[population 2][, etc.]_[allele]"
1. "INVARIANT_![allele]"

If you specify `-v "singletons=1"`, two tags are added:

1. "SINGLETON_all_[allele]" if the allele only exists once across all populations
1. "SINGLETON_[population]_[allele]" if the allele only exists once in the particular population

Arguments:

`singletons`: See above

`exclude`: See above, although this can be a comma-separated list (no whitespace)

`valcols`: Comma-separated list (no whitespace) of column prefixes to use for allele counts, default: AC,AN

Overall, this script allows you to identify shared and population-private
alleles, and filter those alleles for presence or absence in a group.
For instance, you may want to select derived alleles by excluding PanTro,
or to select archaic alleles not found in Africans by excluding Africa and
then looking for "SHARED_[population],[archaic population]_[allele]".

### `extract_by_tag_BED.awk`

This script processes the output of `bcftoolsQuerySiteClassification.awk`
to extract only sites with alleles that match a particular annotation tag.
For instance, if you wanted to extract only sites private to "SouthAsia",
specify `-v "targettag=PRIVATE_SouthAsia"`, which will extract any site
with a "PRIVATE_SouthAsia_[allele]" tag in BED3 format.

Note that you can extract "EXCLUDED_" tags, although if that is not the
target, then any site that is annotated with an "EXCLUDED_" tag is
skipped. "INVARIANT_" tags are also skipped, of course. Also note that
by default, a site where the reference allele has the target tag is
*not* output. Use `-v "incl_ref=1"` to change this.

Arguments:

`targettag`: Tag prefix (except for the allele, of course) to extract

`incl_ref`: Also output sites if the target tag matches the reference allele

### `reformatCVCM_JG.awk`

This script is very custom, and just converts the default output format
of Picard CollectVariantCallingMetrics into a TSV format. The arguments
are strings that are used to parse out some metadata about the variants
input into CVCM.

Arguments:

`prefix`: Prefix of the filename prior to the dbSNP build

`suffix`: Suffix of the filename after the population name

`intermediate1`: String between the dbSNP build and subsampling ID

`intermediate2`: String between the subsampling ID and chromosome number

`intermediate3`: String between the chromosome number and localization of site

This is basically contingent on a filename similar to this:

```[dataset name]_dbSNP[dbSNP build]_VQSRpass_[subsampling ID]_chr[chromosome #]_[localization]_CVCM.variant_calling_summary_metrics```

"localization" here might be "privateEastAsia" which would then be parsed
apart into `sites="private"` and `region="EastAsia"`, representing a VCF of
sites that are private to East Asians. Another example could be
"subsetAfrica" which would be parsed apart into `sites="subset"` and
`region="Africa"`, representing a VCF that just subset out African
samples, didn't filter sites based on population sharing. "subset" and
"private" are the only currently recognized prefixes for localization.

"subsampling ID" is just a label if any subsampling of individuals was
done. For instance, if you wanted to compare allele counts fairly between
populations, you would want to randomly subsample each population down
to the sample number of individuals.

<a name="relatedness" />

## Relatedness

### `extract_families_somalier.awk`

This script recapitulates the trio identification criteria used by
`somalier relate --infer`, but doesn't apply the family-joining algorithm
that has serious trouble with our samples. It also outputs other info
like Parent-Offspring and Full-Sibling pairs (and their relatedness
estimator values) for troubleshooting odd pedigrees and incorrect
branches inferred by --infer.

The inputs are the `[somalier prefix].samples.tsv` file and the
`[somalier prefix].pairs.tsv`, assuming you ran `somalier relate` with
the `SOMALIER_REPORT_ALL_PAIRS` environment variable set to `1` and
with the `--infer` flag.

<a name="sprime" />

## S' (Sprime)


<a name="smcpp" />

## SMC++

### `smcppMergeResults.awk`

This is a simple script to combine the CSVs output by SMC++ together
across runs, labeling based on a parsing of the input file path.
It's not very general, but you could certainly adapt it for your own
filename convention. The point is mainly that the columns are
subset and relabeled so that the R code for plotting is simpler,
and you can co-plot the Ne profiles from each iteration together
(and even co-plot Ne profiles across multiple populations or
parameter sets).
