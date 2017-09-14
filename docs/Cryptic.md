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

