# HumanPopGenScripts
Scripts used for pre- and post-processing of various human population genomics analyses

## Types of analyses:

1. [metadata](#metadata)
1. [admixfrog](#admixfrog)
1. [ADMIXTOOLS input prep](#admixtools)
1. [Archaic allele matching](#archaicmatch)
1. [ArchaicSeeker2.0](#archaicseeker)
1. [De novo assembly](#assembly)
1. [Private allele counting](#privatealleles)
1. [QC](#qc)
1. [Relatedness](#relatedness)
1. [Sprime](#sprime)
1. [SMC++](#smcpp)

<a name="metadata" />

## metadata

### `excludeSamples.awk`

This script takes in a file with sample IDs (one per line) to exclude
(or include) from the second file. The second file must of course
contain a column with sample IDs. The idea is that we can use this
script to filter a list of sample IDs or a metadata file, especially
if we simply maintain a file of sample IDs that are contaminated
or to be pruned due to relatedness, or really anything.

Arguments:

`header`: (optional) Set to 1 if the second file has a header line

`samplecolname`: (optional) Name of the column in the header line corresponding to sample IDs

`negate`: (optional) Set to 1 if you want to *include* the IDs rather than exclude (the default is to exclude)

### `selectSubsamples.awk`

This script takes in a sample metadata file and extracts the sample
IDs with the value(s) indicated in a particular column, possibly
randomly subsampling the samples. This is frequently useful for
extracting the samples belonging to a given population or extracting
equal subsamples from multiple populations.

Arguments:

`idcol`: Column name for the sample ID column in the metadata file

`samplecol`: (optional) Column name to stratify on for subsampling

`subsample`: (optional) Group names to subsample from (comma-separated list)

`sizes`: (optional) Size(s) of subsample(s) to take, must be either one value for all groups or a comma-separated list with one size per group

`seed`: (optional) PRNG seed for reservoir sampling (default: 42)

`selectcol`: (optional) Column name to stratify on for extraction (no subsampling performed with this option)

`select`: (optional) Group names to extract from (comma-separated list)

`subselectsize`: (optional) Size(s) of subsample(s) to take from groups in `select`

If you don't want to perform any subsampling of groups and just extract
all samples in that group, use `selectcol` and `select`. This is the
equivalent of the following on an R data.frame in `tidyverse`:
`df %>% filter([selectcol] %in% c([select])) %>% transmute([idcol])`

If you want to subsample a set of groups, use `samplecol`, `subsample`,
and `sizes` (as well as `seed` if you want). This is similar to the
following on an R data.frame in `tidyverse`:
`df %>% filter([samplecol] %in% c([subsample])) %>% group_by([samplecol]) %>% transmute(id=slice_sample(n=[sizes]))`

It is possible to specify both `select` and `subsample`, in which case
the groups in `select` are extracted and the groups in `subsample` are
subsampled. `select` takes precedence if a group is present in both
lists.

Note that duplicated group names in the lists are not deduplicated, so
the output sample ID list may contain duplicates.

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

<a name="archaicseeker" />

## ArchaicSeeker2.0

### `segToSprimeScore.awk`

This script takes in the `.seg` file output by ArchaicSeeker2 as
well as the output of
`bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'`
on the target/test population VCF to generate a file similar to
the S' `.score` file, except this version generates a tract
per haplotype, not per population.

The idea is that the resulting `.score` file can be filtered
and used to calculate match rates analogous to those of S'
tracts.

The `.seg` file should be sorted and passed as the first
positional argument with a command substitution like this:

```
<(tail -n+2 [.seg file] | \
   sort -k2,2V -k3,3n -k4,4n -k1,1V | \
   cat <(head -n1 [.seg file]) -)
```

Required arguments:

- `chrom`: Comma-separated list of chromosome names

This argument can be used to generate chromosome-specific
`.score` files just like S' produces, or combine them all
into one file, etc.

The output consists of 8 columns (analogous to the S' `.score` file):

1. `CHROM`: Chromosome
2. `POS`: 1-based position of variant (VCF-like)
3. `ID`: Variant ID (VCF-like)
4. `REF`: Reference allele (VCF-like)
5. `ALT`: Alternative allele(s) (VCF-like)
6. `SEGMENT`: String indexing the haplotype, conventionally `[sample ID]_[haplotype index]_[tract index]`
7. `ALLELE`: Allele index of the haplotype at this site (VCF allele index, so non-negative integer)
8. `ASSTATE`: "Best-matching state" of the ArchaicSeeker2.0 tract

Note that `SEGMENT` and `ASSTATE` differ from the canonical S' `.score`
format, as `SEGMENT` is no longer just the positive integer tract index
but now needs to account for the haplotype, and `ASSTATE` replaces the
`SCORE` column, as ArchaicSeeker2.0 does not produce a score, but the
best-matching state may be useful downstream.

Furthermore, note that the results are entirely contingent on the set
of sites output by the `bcftools query` call, so if no site filtering
is done by `bcftools query`, the output will include a line for every
variant in each ArchaicSeeker2.0 tract, regardless of whether it is
informative of archaic ancestry.

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

<a name="qc" />

## QC

### `combine_gtcheck.awk`

This script tidies and combines the output from `bcftools gtcheck -u GT,GT -e 0`
to indicate genotype discordance between all pairs of individuals in
the query and target/"ground-truth" VCFs. The output is sorted in
increasing order of the number of mismatches, so you may want to
re-sort the output in order of the discordance rate (the 5th column)
by piping to `sort -k5,5g`. This generates a rather large file, so
you may want to filter the pairs further for only matching pairs,
but the all pairs results can be useful to detect cryptic matching
pairs, as genotype discordance appears to have a sharp elbow between
matching and non-matching pairs.

### `combine_smplstats.awk`

This script tidies and combines the output from `bcftools +smpl-stats`
across chromosomes. The primary outputs are counts and rates of the
three genotype classes and missingness, though Ti/Tv ratio is also
recalculated and output on a per-sample basis. Thus, the output is
useful for plotting per-sample heterozygosity, missingness, and
Ti/Tv ratio.

### `combine_triostats.awk`

This script tidies and combines the output from `bcftools +trio-stats`
across chromosomes. The primary outputs are various estimates of the
Mendelian error rate, normalized by either the number of valid trio
sites (i.e. no trio member has missing genotype) or the number of
valid trio sites with at least one non-REF allele (i.e. excluding
sites where all trio members are homozygous REF). The calculated
rates include a rate for all putative Mendelian errors, only
homozygous Mendelian errors (i.e. child is homozygous for allele A
and neither parent has allele A), only recurrent Mendelian errors
(i.e. Mendelian errors where the putative error allele is found
in other non-trio samples), and novel singletons (i.e. putative error
alleles found exclusively in the child of a trio).

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

### `archaicMatchSprime.awk`

This script takes a metadata file for the archaics VCF,
the archaics VCF subset to only S' sites, and the S'
`.score` file in order to identify S' sites that match,
mismatch, or are missing from each archaic. This is
the first step toward calculating an archaic match rate
for each S' tract and classifying each tract by it's
origin.

Note that we currently hard-code the expected column
names and archaic group names in the metadata file,
as well as the expected archaic sample IDs. The metadata
file is expected to contain at a minimum two columns
named `Sample` (with the sample ID) and `Region`
(with the archaic hominin species name). These species
names must be `Denisovan` and `Neandertal`. The archaic
sample IDs are expected to be `AltaiNeandertal`,
`Vindija33.19`, `Chagyrskaya-Phalanx`, and `Denisova`.

The output consists of all the S' `.score` file columns plus
8 additional columns:

1. `TractID`: The "conventional" S' tract ID, i.e. `[S' target population]_[chromosome]_[SEGMENT]`
2. `NeandertalMatch`: Categorical variable indicating match state to any Neanderthal
3. `DenisovanMatch`: Categorical variable indicating match state to any Denisovan
4. `NeandertalAlleles`: Comma-separated list of non-missing alleles in Neanderthals
5. `DenisovanAlleles`: Comma-separated list of non-missing alleles in Denisovans
6. `AltaiMatch`: Categorical variable indicating match state to the Altai Neanderthal
7. `VindijaMatch`: Categorical variable indicating match state to the Vindija Neanderthal
8. `ChagyrskayaMatch`: Categorical variable indicating match state to the Chagyrskaya Neanderthal

The categorical variables mentioned here can take on three possible states:

- `match`: The S' allele matches one of the alleles in the query archaic
- `mismatch`: The query archaic genotype is not missing, and the S' allele does not match
- `missing`: The query archaic genotype is missing

The `*Alleles` columns may say `Unk` if all relevant genotypes are missing.

The idea behind the group match (i.e. `NeandertalMatch` for now) is that we
don't want to miss archaic alleles due to sampling limitations when
calculating match rate, so assume that any allele across the three sequenced
Neandertals could be found on an introgressed haplotype. This may inflate
match rates somewhat, but also doesn't deflate them like only using Altai
Neanderthal would.

Required arguments:

- `spop`: S' target population ID (used to construct the conventional S' tract ID)

### `SprimeArchaicMatchRate.awk`

This script takes the output of `archaicMatchSprime.awk`
and summarizes the match indicators into match rates and
counts of ascertainable ("Good") sites.

For now, the groups for which match rates and "Good" site
counts are determined are hard-coded in the `BEGIN` block
as being all Neanderthals (`Neandertal`), the Altai Denisovan
(`Denisovan`), and each of the Neanderthals individually
(`Altai`, `Vindija`, and `Chagyrskaya`).

The output consists of 6 core columns, and 2 columns per
group, so 16 columns as of the writing of this README:

1. `CHROM`: Chromosome
2. `START`: 1-based S' tract start position
3. `END`: 1-based S' tract end position
4. `TractID`: S' tract ID
5. `SNPLEN`: Number of S' variants in the tract
6. `SCORE`: S' score of the tract

The group columns are:

7. `[NDAVC]matchrate`: Proportion of matching sites to ascertainable sites for group [NDAVC]
8. `[NDAVC]good`: Number of ascertainable sites for group [NDAVC]

This pair of columns is repeated for each value in `[NDAVC]`,
representing the aforementioned 5 hard-coded group.s

### `extract_Sprime_arcmatch_sites.awk`

This script takes in the S' match rates file produced by
`SprimeArchaicMatchRate.awk` as well as the S' archaic matches
file produced by `archaicMatchSprime.awk` to identify a set
of sites from each S' tract, optionally only those that
match a given archaic hominin.

This step is analogous to the first tag SNP filtering step in
the adaptive introgression haplotype identification pipeline
of [Gittelman et al. 2016 Current Biology](https://doi.org/10.1016/j.cub.2016.10.041)
where only tag SNPs matching the Altai Neanderthal are selected
for evaluation of r^2.

The output of this script can be used to select SNPs within the
target population (or superpopulation) for calculation of r^2.
It is simply a 2-column TSV of `CHROM` and `POS`, as expected by
the `-R` and `-T` arguments of the `bcftools` commands.

Required arguments:

- `source`: Which archaic origin to select tracts from (e.g. Neandertal or Denisovan)

Optional arguments:

- `criteria`: Which match rate criteria to use (i.e. Browning or PFR for now)
- `only_matches`: Only output sites in source tracts that also match the source

Criteria definitions:

`Browning` (i.e. [Browning et al. 2018 Cell](https://doi.org/10.1016/j.cell.2018.02.031)):

- At least 30 ascertainable sites in both Altai Neanderthal and Altai Denisovan
- `Neandertal` if Altai Neanderthal match rate >= 0.6 and Altai Denisovan match rate <= 0.4
- `Denisovan` if Altai Neanderthal match rate <= 0.3 and Altai Denisovan match rate >= 0.4

`PFR` (my own ad hoc criteria):

- At least 30 ascertainable sites in all 3 Neanderthals and the Altai Denisovan
- `Neandertal` if combined Neanderthal match rate >= 0.3 and Altai Denisovan match rate <= 0.3
- `Denisovan` if combined Neanderthal match rate <= 0.3 and Altai Denisovan match rate >= 0.3
- `Ambiguous` if combined Neanderthal match rate >= 0.3 and Altai Denisovan match rate >= 0.3

### `SprimeArchaicAF.awk`

This script takes the S' `.score` file and the output of
`bcftools query -H -f '%CHROM:%POS[\t%GT]\n'` on a VCF
containing samples from a query population to calculate
the frequency of the putative S' archaic allele at each
site. Downstream, this can be processed into an estimate
of the S' tract frequency in the query population.

Note that you'll want to run `bcftools query` on the
subset of sites within the S' `.score` file.

The output consists of 7-8 columns:

1. `CHROM`: Chromosome
2. `POS`: 1-based position of the variant
3. `TractID`: S' tract ID, typically taking the form `[S' target population]_[chromosome]_[SEGMENT]`
4. `ArchaicAlleleCount`: Number of query population alleles matching the S' putative archaic allele
5. `TotalAlleleCount`: Number of non-missing alleles in the query population
6. `ArchaicAlleleFrequency`: The quotient of the previous two columns
7. `Population`: ID of the query population

An additional column may be inserted before the `ArchaicAlleleCount`
column if the `allele` argument is set:

- `ArchaicAllele`: The S' allele at the current site (not the allele index)

This is mainly useful for diagnostics, as it breaks compatibility with
downstream tools like `SprimeTractMedianAF.awk`.

Required arguments:

- `pop`: ID of the query population
- `spop`: ID of the S' target population (this is used to construct the conventional S' TractID)

Optional arguments:

- `header`: Flag to indicate whether to output a header (default: no header)
- `allele`: Flag to indicate whether to output the additional `ArchaicAllele` column as the fourth column (default: no additional column)

### `SprimeTractMedianAF.awk`

This script takes the output of `SprimeArchaicAF.awk` and
estimates the frequency of the entire tract as the median
of the site-wise frequencies along that tract. The minimum,
maximum, and number of sites considered are also output
as diagnostics for the breakdown of the tract haplotype
across individuals. For instance, if min, median, and max
are all similar then the tract haplotype hasn't broken
down much, whereas if min is much smaller than median and
max, the tract haplotype has likely broken down toward the
edges, and if min and median are both smaller than max, then
the tract haplotype is likely very fragmented. The median
is a poor estimate when the number of sites is small, so
beware those cases as well.

The output consists of 9 columns:

1. Chromosome
2. Tract start position (0-based, BED-like)
3. Tract end position (1-based, BED-like)
4. TractID
5. Query population (i.e. the population you want tract frequencies for)
6. Median S' allele frequency
7. Minimum S' allele frequency
8. Maximum S' allele frequency
9. Number of sites considered

### `SprimePerSampleTracts.awk`

This script takes the S' `.score` file and the output of
`bcftools query -H -f '%CHROM:%POS[\t%GT]\n'` on the VCF of modern
samples, and outputs the state of each sample at each site
in terms of genotype match to the S' putative archaic allele.

The idea is that you can run `bcftools query` on the subset of
S' sites, and this script then projects the genotypes of each
individual onto the S' tracts, polarizing the genotypes against
the S' allele state.

The output consists of four columns:

1. `CHROM:POS`
2. `Individual`: The sample ID of the individual
3. `TractID`: The ID of the S' tract
4. `TractState`: The polarized genotypic state at that site

`Individual` may be replaced by `Haplotype` if the `phased` flag is set,
in which case it represents a haplotype ID.

`TractState` can take on up to three values:

- `homozygous` (i.e. homozygous for the S' putative archaic allele)
- `heterozygous` (i.e. heterozygous for the S' putative archaic allele)
- `homozygous_nonarchaic` (i.e. homozygous for the putative modern allele)

Arguments:

- `allout`: A flag to include `homozygous_nonarchaic` sites in the output (default: only output `homozygous` and `heterozygous` sites)
- `header`: A flag to output a header line (default: no header)
- `phased`: A flag to indicate that the input is phased, so haplotypes should be analyzed instead of genotypes (default: analyze genotypes)

### `SprimeTractBED.awk`

This script takes the output of `SprimePerSampleTracts.awk`, which is
a table of positions in S' tracts and their genotypic states in each
individual, and outputs the boundaries of these tracts in each individual.
In the genotypic case, a tract can be extended by a homozygous or
heterozygous variant, so it may represent a composite of two haplotypes.

The output is a BED6 file with semicolon-delimited tag-value pairs in the
name column (column #6):

- `TractID`: S' tract ID, formed as `[population]_[chromosome]_[SEGMENT]`
- `Individual`: Sample/individual ID
- `State`: Polarized genotypic state along the tract (one of the values in column 4 of the input, or "mixed" if that state changes within the tract)
- `HomSprimeSites`: Count of sites with the "homozygous" state (i.e. homozygous for the putative archaic allele)
- `HetSprimeSites`: Count of sites with the "heterozygous" state (i.e. heterozygous for the putative archaic allele)

### `SprimeTractBEDtoLengths.awk`

This is a very basic script that takes a file of sample IDs and
the output of `SprimeTractBED.awk` to calculate an estimate of
total archaic sequence per individual. Output columns are:

1. `Chromosome`
2. `Population`
3. `TractID` (see above)
4. `Sample` (i.e. Individual above, or Haplotype if `phased` is set)
5. `TRACTLEN` (length of the tract)
6. `ARCMOD` (number of heterozygous sites in the tract, i.e. HetSprimeSites)
7. `ARCARC` (number of homozygous archaic sites in the tract, i.e. HomSprimeSites)

Required arguments:

- `pop`: Used as a prefix for the output TractID

Optional arguments:

- `header`: Output a header line
- `phased`: Consider input as phased, and generate haplotypic tracts

### `SprimeTractGeneList.awk`

This script reformats the output of `bedtools intersect -wao`
between a BED6 of S' tracts produced by `SprimeTractBED.awk`
(fed to `bedtools intersect` as the `-a` argument) and a GFF3
to indicate the set of features overlapping each S' tract.
If multiple features overlap an S' tract, they are concatenated
into a comma-separated list such that each S' tract has a single
output record.

Arguments:

- `feature`: The GFF3 feature type (i.e. column 3) to check for overlaps (default: gene)
- `tag`: The tag of the GFF3 feature to extract a name from (default: ID)
- `trim`: Flag to trigger trimming of "[.][0-9]+(_[0-9]+)?$" from the end of the tag value

So if you want the canonical gene names in the list rather than e.g.
Ensembl gene IDs, set `tag` to `gene_name`.

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
