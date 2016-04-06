##Discovering Ancestry
1000Genomes contained data from European, African, East Asian and American samples. The American samples
are expected to be admixed with a lot of European ancestry. We therefore expect principle component analysis to give us
3 clusters, corresponding to European, East Asian and African with the American samples either in or very close to the European 
cluster. The command to run is
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/1000G.snps.nochr.vcf.gz > test_pca
```
or if you want to use unrelated samples use the subsetting option
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/1000G.snps.nochr.vcf.gz -S phase1.unrel > test_pca.unrel
```
Alternatively we can project samples onto 1000Genomes phase3 princple components using
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -w data/1000G.snps.nochr.vcf.gz > test_pcaproj
```

Plotting the projections of all samples onto the first three principle components gives

![alt text](https://github.com/Illumina/akt/blob/master/docs/test_pca12.png)
![alt text](https://github.com/Illumina/akt/blob/master/docs/test_pca32.png)

The colours are assigned based on the known ancestry provided by 1000Genomes, showing that the first principle components
do a good job of classifying samples based on ancestry for African and East Asian samples. European samples
are also tighly clustered with American samples close by. This is indicative of an admixed population.

If we didn't know the ancestry we could use a clustering algorithm to classify the samples into different sub-populations.
```
./akt cluster test_pca -c 2-3 -a 2 -d 1 --density-plot > test_dplot
./akt cluster test_pca -c 2-3 -a 2 -d 1 -p -1 -D 700 > test_allclusters
```
The `-a 2` option uses a density based clustering method. 
The first command locates density peaks, which are cluster centres, the `-d 1` counts neighbours in a ball of
radius 1 around each point to estimate local density. The second command does clustering with parameters
deduced from `phase1_dplot`. We use a cutoff for peaks of 700 and set the minimum density to -1, so that
all points are classified.

![alt text](https://github.com/Illumina/akt/blob/master/docs/test_cluster.png)

We classify the data into 5 groups corresponding the the 1000G superpopulations. If we don't insist on
classifying all the data we can localise the clusters to where the data is densest.

```
./akt cluster test_pca -c 2-3 -a 2 -d 1 -p 0 -D 700 > test_allclusters1
```

![alt text](https://github.com/Illumina/akt/blob/master/docs/test_cluster1.png)

Cluster0 is the group of samples which were unclassified.
