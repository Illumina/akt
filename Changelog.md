##Changelog

##2016.06.03
* kin re-factoring
 * added `-M 0/1/2` for plink/kin/grm kinship calculations
 * frequencies now read from `-F` argument and `-R/-T` now behave as bcftools
 * moved kinship calculations into their own class
 * reduced memory footprint (at the cost of non-deterministic output ordering)

##2016.06.03
* kin command defaults to calculating allele frequency from data
* changed allele frequency vcf argument to -F

##2016.05.06
* added a check for 0 samples in akt pca

##2016.04.08
* initial public release

