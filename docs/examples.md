# Example Workflow
## Test Data

First we need some data. To see some interesting behaviour we should have a lot of samples and some of them should be related. 
We use publicly available [1000 Genomes Project](http://www.1000genomes.org/) data. We will work with a BCF containing 433 high-coverage samples with a mix of different ethnicities as well as 129 mother-father-child and 9 parent-child duos. See [this script](https://gist.github.com/jaredo/4206a09eedc7a0fed3f09ca756af0919) for a description of how the BCF was generated. 

We can download the data via:
```
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi
```
AKT only works with autosomal DNA i.e. not X and Y, so you should not use these chromosomes.
## Discovering Cryptic Relations
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

## Discovering Ancestry
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
