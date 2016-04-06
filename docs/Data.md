##Test Data

First we need some data. To see some interesting behaviour we should have a lot of samples and some of them should be related. 
We use publicly available [1000 Genomes Project](http://www.1000genomes.org/) data. We will work with a BCF containing 433 high-coverage samples with a mix of different ethnicities as well as 129 mother-father-child and 9 parent-child duos. See [this script](https://gist.github.com/jaredo/4206a09eedc7a0fed3f09ca756af0919) for a description of how the BCF was generated. 

We can download the data via:
```
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi
```
AKT only works with autosomal DNA i.e. not X and Y, so you should not use these chromosomes.
