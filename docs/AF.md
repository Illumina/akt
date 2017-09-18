## Allele Frequency and Correlation

We can use the PCA analysis above to calculate population specific allele frequencies. We choose
the isolated cluster 1 (East Asian). We can either subsample the input bcf or use the subsetting option `-S`. To
save time we only do chromosome 20:
```
grep Cluster1 test_allclusters | awk '{print $4}' > test.Cluster1
./akt stats ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -a "C1" -O b -o ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.sitesx.bcf -S test.Cluster1 -r 20 -x
```
the `-a` option will label the allele frequency "C1_AF" in the output instead of "AF".
To calculate LD score we leave out the -x option and look at say, 10K flanking bases, this will take much longer, so
try just a short chromosome for this test.
```
./akt stats ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -a "EAS" -O b -o sites.EAS.bcf -S test.Cluster1 -r 20 -F 10000
bcftools index sites.EAS.bcf
``` 
