##kin
```
replace_kin_run

```
* -k : Output only pairs with kinship coefficient more than this.
* -u : Don't normalize IBD scores
* -c : Calculate allele frequency from data
* thin_option
* targetfile_option
* nthread_option
* aftag_option
* regions_option
* regionsfile_option
* samples_option
* samplesfile_option
* pairfile_option
* maf_option


Run the kinship calculation by giving akt a multi-sample vcf/bcf file:
```
~/ancestry_tools$./akt kin multisample.bcf -R data/1000G.snps.nochr.vcf.gz -n 32 > allibd
```
the output `allibd` is 7 columns in the format
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

