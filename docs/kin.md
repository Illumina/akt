##kin

Calculates kinship coefficients (and other related metrics) from multi-sample VCF/BCFs. Can be used to detect (closely) related or duplicated samples.

**-k** *value* Only output pairs with kinship coefficient greater than *value*  
**-F** *FILE* a file containing population allele frequencies to use in kinship calculation  
**-M** *0/1/2* type of estimator.  0:[plink (default)](https://www.cog-genomics.org/plink2/ibd) 1:[king-robust](http://bioinformatics.oxfordjournals.org/content/26/22/2867.full) 2:[genetic-relationship-matrix](http://cnsgenomics.com/software/gcta/estimate_grm.html)  
aftag_option
nthread_option
regionsfile_option
regions_option
targetfile_option
target_option
samplesfile_option
samples_option


Run the kinship calculation by giving akt a multi-sample vcf/bcf file:

```
$ akt kin multisample.bcf -R data/wgs.grch37.vcf.gz -n 32 > kin.txt
```
###Choice of estimator


The default algorithm (`-M 0`) used to calculate IBD is taken from [PLINK](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/) with some minor changes. It outputs the following seven column format:

```
ID1 ID2 IBD0 IBD1 IBD2 KINSHIP NSNP
```

As with PLINK, we set KINSHIP = 0.5 * IBD2 + 0.25 * IBD1. Our IBD values may slighly differ to PLINK's (by desing) due to the following differences:

* No 'bias correction' since allele frequencies are assumed to be accurate
* Normalization as follows:
	* if IBD0 > 1: IBD0 = 1, IBD1 = 0, IBD2 = 0
	* if IBD1 < 0: IBD1 = 0
	* if IBD2 < 0: IBD2 = 0
	* norm = IBD0 + IBD1 + IBD2
	* IBD0 /= norm, IBD1 /= norm; IBD2 /= norm;
* We do **not** follow PLINK which forces IBD to obey consistency conditions - this affects the clustering that is required for the `relatives` code.

The second method (`-M 1`) uses the robust kinship coefficent estimate describing in the [KING paper](http://bioinformatics.oxfordjournals.org/content/26/22/2867.full). This may be preferable when your cohort has large amounts of population structure. The IBD estimates and output format are the same as for `-M 0`.

The third method (`-M 2`) outputs the genetic relationship matrix used by [GCTA](http://cnsgenomics.com/software/gcta/) (among other software).  The output format is `ID1 ID2 GR` where `GR` is the genetic correlation between samples `ID1` and `ID2`. This matrix is important for fitting linear mixed-effect models but we have found methods 0 and 1 more appropriate for sample QC purposes.

