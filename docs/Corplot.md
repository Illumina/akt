##Correlation Matrix
We can visually identify linkage disequilibrium by plotting a correlation matrix.

```
./akt LDplot ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20:0-200000 > test_sigma
```
We only look at a 200kb region of chromosome 20 since calculaing very large correlation matrices is very time 
and memory consuming.
```
plot 'test_sigma' matrix with image
```

![alt text](https://git.illumina.com/Bioinformatics/akt/blob/master/docs/test_sigma.png)

A value near 1 indicates very strong positive correlation (variants always present together)
and a value near -1 indicates string negative correlation (variants never present together). 
Looking at the plot above showing variants 800 to 1200 in the region 20:0-200000 and their correlation
we see obvious LD blocks.
