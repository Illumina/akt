##kin
```
replace_kin_run

```
* -k : Only output pairs with kinship coefficient greater than this.
* -F : a file containing population allele frequencies to use in kinship calculation
* -M : type of estimator. 0:[plink (default)](https://www.cog-genomics.org/plink2/ibd) 1:[king-robust](http://bioinformatics.oxfordjournals.org/content/26/22/2867.full) 2:[genetic-relationship-matrix](http://cnsgenomics.com/software/gcta/estimate_grm.html)
* aftag_option
* nthread_option
* regionsfile_option
* regions_option
* targetfile_option
* target_option
* maf_option
* thin_option
* samples_option
* samplesfile_option


Run the kinship calculation by giving akt a multi-sample vcf/bcf file:

```
$ akt kin multisample.bcf -R data/wgs.grch37.vcf.gz -n 32 > kin.txt
```

the output `kin.txt` is 7 columns in the format

```
ID1 ID2 IBD0 IBD1 IBD2 KINSHIP NSNP
```

The algorithm used to calculate IBD is taken from [PLINK](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/) with some minor changes:

* No 'bias correction' since allele frequencies are assumed to be accurate
* Normalization as follows:
	* if IBD0 > 1: IBD0 = 1, IBD1 = 0, IBD2 = 0
	* if IBD1 < 0: IBD1 = 0
	* if IBD2 < 0: IBD2 = 0
	* norm = IBD0 + IBD1 + IBD2
	* IBD0 /= norm, IBD1 /= norm; IBD2 /= norm;
* We do NOT follow PLINK which forces IBD to obey consistency conditions - this affects the clustering that is required for the `relatives` code.
* KINSHIP = 0.5 * IBD2 + 0.25 * IBD1;

