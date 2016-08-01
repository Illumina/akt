#Example Workflow
##Test Data

First we need some data. To see some interesting behaviour we should have a lot of samples and some of them should be related. 
We use publicly available [1000 Genomes Project](http://www.1000genomes.org/) data. We will work with a BCF containing 433 high-coverage samples with a mix of different ethnicities as well as 129 mother-father-child and 9 parent-child duos. See [this script](https://gist.github.com/jaredo/4206a09eedc7a0fed3f09ca756af0919) for a description of how the BCF was generated. 

We can download the data via:
```
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi
```
AKT only works with autosomal DNA i.e. not X and Y, so you should not use these chromosomes.
##Discovering Cryptic Relations
First we attempt to identify cryptic relations. Run
```
./akt kin ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/wgs.grch37.vcf.gz -n 4 > test_kin
```
To calculate average ibd for all pairs in the dataset. Use more than 4 processors if you have them! `test_kin` contains 
IBD and kinship coefficients for all 93528 possible pairs, most of which are not closely related.

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_kin.png)


The `relatives` tool identifies hidden relatives.
```
./akt relatives test_kin -p test -g > test_relatives
```
First have a look at the graph of all relatives (the graph is created using [GraphViz](http://www.graphviz.org/))
```
neato test.allgraph -Tpng -O
eog test.allgraph.png
```

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test.allgraph.png)

`relatives` also attempts to provide more information about the pedigree structure in each relative group that it finds. 
These are output as graph files `test.Fam*.graph`. Sometimes there is not enough information to tell 
which sample is the parent and which is the offspring. In this case an arbitrary choice is made in the .fam file and the graph file 
represents this as a double arrow.

```
dot test.Fam0.graph -Tpng -O
```
![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test.Fam0.graph.png)

some families can be resolved correctly
```
dot test.Fam134.graph -Tpng -O
```
![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test.Fam134.graph.png)

Running the command
```
grep '^Type' test_relatives
```
gives a list of sample pairs together with their family id and the best fit to the degree of relationship between the two 
samples. This is in good accord with the list of relatives found by 1000Genomes:
``` 
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt 
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/README.sample_cryptic_relations.
```

We can pull out an unrelated subset of the 1000G phase1 individuals using
```
grep '^Unrel' test_relatives | awk '{print $2}' > test.unrel
```
and use the `-S` argument as we proceed to avoid creating large temporary files. These samples
should be at worst 3rd order relatives.

This dataset contains duplicate samples which can be found with
```
grep '^Dup' test_relatives
```
 
Duplicates, parent/child relations and sibling relationships should be quite well determined by `akt relatives` but 
higher orders may be difficult to distinguish e.g. 2nd order from 3rd order relations and especially 3rd order relations from
higher relations. `akt relatives` uses a very simple approach - assign the relationship to the category it is closest
to (Euclidean distance) in (IBD0,IBD1) space. More complicated approaches may be implemented in future - e.g. mixture
distributions, but for now expect nuclear families to
be well resolved but a network of cousins is probably going to be missed.

##Profiling Mendelian inheritance
The `relatives` analysis in the previous section provided us with a set of likely pedigrees:
```
$ head test.fam 
Fam0	NA18909	0	0	0	0
Fam0	NA18911	NA18909	0	0	1
Fam1	HG00533	0	0	0	0
Fam1	HG00534	0	0	0	0
Fam1	HG00535	HG00533	HG00534	0	2
Fam2	HG00536	0	0	0	0
Fam2	HG00537	0	0	0	0
Fam2	HG00538	HG00536	HG00537	0	2
Fam3	HG00619	0	0	0	0
Fam3	HG00620	0	0	0	0
```
We can profile Mendelian inheritance in these pedigrees via the `akt mendel` command. This allows us to double-check the discovered pedigrees are correct ie. they have very few Mendelian inconsistent genotypes. It also provides a useful quality-control metric, in that variants with high amounts of Mendel inconsistency are likely to be problematic.

First we run the `mendel` subcommand using the `.fam` file and the genotype BCF as input:
```
$ akt mendel -p test.fam ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf > mendel.out
Read 405 individuals from n433.fam
Found 405 in both the pedigree and bcf.
Found 129 trios and 9 duos
Reading input from ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
34258858 variants checked
```

The output is very granular, with one line per every possible parental genotype configuration. We can summarise the output using a provided R script:

```
$ Rscript  akt_dir/scripts/mendel.R mendel.out

Trio summary:
 DAD MUM         RR       RA       AA error_rate     het_rate
  AA  AA        136     4429 79364648 0.00575160  0.005580249
  AA  RA      12119 19194531 18892662 0.03180897 50.380256210
  AA  RR      32000 18995403    14511 0.24425591 99.755744092
  RA  AA      11356 19123054 19040795 0.02974706 50.092865251
  RA  RA   18954076 37544501 18483645 0.00000000 50.071203545
  RA  RR   56846315 54590875    20426 0.01832625 48.979044196
  RR  AA      29879 19097905    15895 0.23910764 99.760892355
  RR  RA   56791975 54948701    21004 0.01879356 49.165958314
  RR  RR 3823192009  1512923     3671 0.03965254  0.039556556

Duo summary:
 DAD MUM        RR       RA      AA  error_rate het_rate
  AA   .      4628  2906054 6726522 0.048022227 30.15453
  RA   .   7237430 10001156 2937120 0.000000000 49.57029
  RR   . 264510589  7248782    4355 0.001602495  2.66731
```
Everything looks consistent here, there are low (<0.5%) rates of Mendel inconsistencies and the transmission rate of heterozygous variants is close to 50% where appropriate. This is not suprising since this data set is of high quality.
##Discovering Ancestry
1000Genomes contained data from European, African, East Asian and American samples. The American samples
are expected to be admixed with a lot of European ancestry. We therefore expect principle component analysis to give us
3 clusters, corresponding to European, East Asian and African with the American samples either in or very close to the European 
cluster. The command to run is
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/wgs.grch37.vcf.gz > test_pca
```
or if you want to use unrelated samples use the subsetting option
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/wgs.grch37.vcf.gz -S phase1.unrel > test_pca.unrel
```
Alternatively we can project samples onto 1000Genomes phase3 princple components using
```
./akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -W data/wgs.grch37.vcf.gz > test_pcaproj
```

Plotting the projections of all samples onto the first three principle components gives

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_pca12.png)
![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_pca32.png)

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

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_cluster.png)

We classify the data into 5 groups corresponding the the 1000G superpopulations. If we don't insist on
classifying all the data we can localise the clusters to where the data is densest.

```
./akt cluster test_pca -c 2-3 -a 2 -d 1 -p 0 -D 700 > test_allclusters1
```

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_cluster1.png)

Cluster0 is the group of samples which were unclassified.
##Allele Frequency and Correlation

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

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_sigma.png)

A value near 1 indicates very strong positive correlation (variants always present together)
and a value near -1 indicates string negative correlation (variants never present together). 
Looking at the plot above showing variants 800 to 1200 in the region 20:0-200000 and their correlation
we see obvious LD blocks.
##Calculating Admixture Fractions
It is possible to linearly transform principle components to a new basis in which the co-ordinates are the 
admixture fractions of different ancestral populations: see 
the paper of [Zheng and Weir](http://www.sciencedirect.com/science/article/pii/S0040580915000891). We provide a very
simple method to find this linear transformation. For example by looking at the plot below (projecting onto the first two
1000G principle components and removing the SAS samples for clarity.)

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_pcaproj12.png)

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

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_admix.png)

The transformation specified in
data/1000G.pca_to_admix works reasonably well for admixture in 1000G superpopulations
(AFR, AMR, EAS, EUR, SAS) if (and ONLY if) the PCA is done by projection using the
`-w` option and the file `wgs.grch37.vcf.gz` or `wgs.hg19.vcf.gz`
```
./akt admix test_pcaproj -c 2-6 -C data/1000G.pca_to_admix > test_alladmix
```

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_alladmix.png)


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
![alt text](https://github.com/Illumina/akt/blob/master/docs/test_metafreq.png)

The QF score is correlated with allele frequency difference as it should be. Of course, the QF score
is much more powerful for significance testing. As well as differentiating
sites in unrelated populations this tool could also be of use for quality control of different bioinformatics
pipelines. Systematic errors made in one but not the other would be highly statistically significant and so have 
large QF and QX values.
