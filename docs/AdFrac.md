##Calculating Admixture Fractions
It is possible to linearly transform principle components to a new basis in which the co-ordinates are the 
admixture fractions of different ancestral populations: see 
the paper of [Zheng and Weir](http://www.sciencedirect.com/science/article/pii/S0040580915000891). We provide a very
simple method to find this linear transformation. For example by looking at the plot below (projecting onto the first two
1000G principle components and removing the SAS samples for clarity.)

![alt text](https://git.illumina.com/Bioinformatics/akt/blob/master/docs/test_pcaproj12.png)

We can see that the centre of the African cluster is near (-49,-4), the European cluster near (-13,27) and the East Asian
cluster near (-25-30). Let these points represent 100% African, European and East Asian respectively. To transform
the pca plot into a plot of admixture fractions - say %African and %European form the following table
```
cat centre.txt
40 -4 AFR 1 0
-12 21 EUR 0 1 
-19 -27 EAS 0 0
```
then run 
```
./akt admix test_pcaproj -c 2-3 -C centre.txt > test_admix
```

![alt text](https://git.illumina.com/Bioinformatics/akt/blob/master/docs/test_admix.png)

The transformation specified in
data/1000G.pca_to_admix works reasonably well for admixture in 1000G superpopulations
(AFR, AMR, EAS, EUR, SAS) if (and ONLY if) the PCA is done by projection using the
`-w` option and the file `1000G.snps.nochr.vcf.gz` or `1000G.snps.withchr.vcf.gz`
```
./akt admix test_pcaproj -c 2-6 -C data/1000G.pca_to_admix > test_alladmix
```

![alt text](https://git.illumina.com/Bioinformatics/akt/blob/master/docs/test_alladmix.png)


