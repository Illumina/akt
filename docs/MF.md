##Checking Site Frequency Distributions
Given two populations, or even the same population called using two different pipelines, we can
ask if the observed allele frequencies are consistent with each other. For example: say in cohort 1
at a particular site we observed 20 variants out of 200 possible alleles while in cohort 2 we observed 15 out of 100.
Is the difference in allele frequency 0.1 versus 0.15, significant given the cohort sizes? 
We provide two statistical tests: Fisher's exact test and the chi-squared test.

For example - we can try to find at which sites European and East Asian samples have statistically different allele frequencies.
First we subset the data
```
bcftools view ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -S EAS.samples -r 20 -O b -o EAS.20.bcf --force-samples
bcftools view ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -S EUR.samples -r 20 -O b -o EUR.20.bcf --force-samples
bcftools index EAS.20.bcf
bcftools index EUR.20.bcf
```
Then run `akt metafreq`
```
./akt metafreq EAS.20.bcf EUR.20.bcf -O b -o mf.bcf
```

The file `mf.bcf` is a site only VCF containing -log10(p-value) for each site where p-values
are calculated using Fisher's exact (QF) and chi-squared (QX) tests. Higher values indicate
more significant differences. A good predictor
of population differences is simply the difference in allele frequency.
Observing AF=0.99 in one group and AF=0.05 in another is a sure sign those groups differ at that site.
![alt text](https://raw.githubusercontent.com/Illumina/akt/metafreq/docs/test_metafreq.png)

The QF score is correlated with allele frequency difference as it should be. Of course, the QF score
is much more powerful for significance testing. As well as differentiating
sites in unrelated populations this tool could also be of use for quality control of different bioinformatics
pipelines. Systematic errors made in one but not the other would be highly statistically significant and so have 
large QF and QX values.
