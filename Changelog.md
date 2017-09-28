## Changelog

## 2017.08.28
* pedphase now handles the PS tag correctly
* better error handling when reading the pedigree
* more informative help for pedphase

## 2017.08.19
* added VCF4.3 pedigree parsing

## 2017.08.09
* added unrelated command
* relatives command no longer prints unrelated (see above)
* added the pedphase command

## 2017.05.02
* added --assume-homref flag to pca

## 2017.03.01
* mendel.R tweak
* stdout tweaks	
* better error handling in pedigree (.fam) reader

## 2016.08.01
* removed akt ibd

## 2016.06.03
* kin re-factoring
 * added `-M 0/1/2` for plink/kin/grm kinship calculations
 * frequencies now read from `-F` argument and `-R/-T` now behave as bcftools
 * moved kinship calculations into their own class
 * reduced memory footprint (at the cost of non-deterministic output ordering)
 * removed pairfile option
 * users must provide one of -F/-R/-T or use --force

## 2016.06.03
* kin command defaults to calculating allele frequency from data
* changed allele frequency vcf argument to -F

## 2016.05.06
* added a check for 0 samples in akt pca

## 2016.04.08
* initial public release

