# akt - ancestry and kinship toolkit
##License

AKT is freely available under the [GPL3 license](https://github.com/Illumina/agg/blob/master/LICENSE). AKT is not commercially supported.

AKT relies on HTSlib and Eigen. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a header-only library for matrix algebra released under the [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/). [HTSlib](http://www.htslib.org/) is a library for efficently parsing vcf/bcf files released under the [MIT/Expat License](http://choosealicense.com/licenses/mit/).

##Installation instructions

```
git clone git@git.illumina.com:Bioinformatics/akt.git
make
```

##Using akt
akt uses the syntax
```
./akt COMMAND <options>
```
To see a list of available tools use
```

Program:	akt (Ancestry and Kinship Tools)
Version:	ba2655b
Contact:	rarthur@illumina.com

Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.

Usage:	akt <command> [options]

	pca                      principal component analysis
	kin                      detect average IBD sharing
	relatives                discover pedigrees
	ibd                      detect segments shared IBD
	mendel                   profile Mendelian inhertiance and inconsistencies in known pedigrees
	cluster                  perform cluster analyses
	LDplot                   output correlation matrix
	stats                    calculate AF and LD metrics

```
##kin
```
Calculate IBD stats from a VCF
Usage:
./akt kin in.bcf -T sites.vcf.gz
Expects input.bcf to contain genotypes.
User must specify a file containing allele frequencies with -R!
	 -k --minkin:			threshold for relatedness output (none)
	 -u --unnorm:			If present don't normalize
	 -c --calc:			calculate AF from data
	 -h --thin:			keep every t variants
	 -T --targets-file:		intersecting VCF
	 -n --nthreads: 		num threads
	 -a --aftag:			allele frequency tag
	 -r --regions:			chromosome region
	 -R --regions-file:		list of regions, file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
	 -f --pairfile:			file containing sample pairs
	 -m --maf:			minimum MAF

```
* -k : Output only pairs with kinship coefficient more than this.
* -u : Don't normalize IBD scores
* -c : Calculate allele frequency from data
* -h : Hop this many markers forward from last accepted marker. 
* -T : Indexed VCF file containing intersecting sites and relevant site info. 
* -n : Number of threads to use. 
* -a : Allele frequency tag e.g. "TAG" reads or writes "TAG_AF" in the target VCF. 
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -s : Comma-separated list of samples to include or exclude if prefixed with "^".
* -S : File of sample names to include or exclude if prefixed with "^". One sample per line.
* -f : File containing pairs of samples to compute, 2 column: Sample1 Sample2. 
* -m : Minimum MAF a site to be counted. When reading: if AF > 0.5 checks 1-AF > m. 


Run the kinship calculation by giving akt a multi-sample vcf/bcf file:
```
~/ancestry_tools$./akt kin multisample.bcf -R data/1000G.snps.nochr.vcf.gz -n 32 > allibd
```
the output `allibd` is 7 columns in the format
```
ID1 ID2 IBD0 IBD1 IBD2 KINSHIP NSNP
```
The algorithm used to calculate IBD is taken from [PLINK](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/) with some minor changes:
* No 'bias correction' since allele frequencies are assumed to be accurate
* Normalization as follows:
	* if IBD0 > 1: IBD0 = 1, IBD1 = 0, IBD2 = 0
	* if IBD1 < 0: IBD1 = 0
	* if IBD2 < 0: IBD2 = 0
	* norm = IBD0 + IBD1 + IBD2
	* IBD0 /= norm, IBD1 /= norm; IBD2 /= norm;
* We do NOT follow PLINK which forces IBD to obey consistency conditions - this affects the clustering that is required for the `relatives` code.
* KINSHIP = 0.5 * IBD2 + 0.25 * IBD1;

##relatives
```
Discover relatives from IBD
Usage:
./akt relatives ibdfile
	 -k --kmin:			threshold for relatedness (0.05)
	 -i --its:			number of iterations to find unrelated (10)
	 -g --graphout:			if present output pedigree graph files
	 -p --prefix:			output file prefix (out)
arrow types     : solid black	= parent-child
                : dotted black	= siblings
                : blue 		= second order
                : red		= duplicates
                : directed	= from parent to child

```
* -k : Only keep links with kinship above this threshold (searches in this set for duplicate, parent-child and sibling links).
* -i : Iteration parameter for unrelated set output.
* -p : Prefix for output files.
* -g : If present output graphviz files. These can be visulaised using e.g. `neato -Tpng -O out.allgraph` or for family pedigrees `dot -Tpng -O out.Fam0.graph`.
```
./akt relatives allibd -g > allrelatives
```
The output contains duplicates, families, relationship types and unrelated individuals
```
grep ^Dup allrelatives
Dup0 Sample0
Dup0 Sample1
...
grep ^Fam allrelatives
Fam0 Sample2
Fam0 Sample3
...
...
grep ^Type allrelatives
Type Fam0 Sample2 Sample3 Parent/Child
...
grep ^Unrel allrelatives
Sample0
Sample2
...
```
The file `out.allgraph` can be viewed with gviz e.g. `fdp out.allgraph -Tpng -O` and the families can be viewed using 
e.g. `dot out.Fam0.graph -Tpng -O`. The parent child relations are also recorded in PLINK fam format in `out.fam`. If
e.g. a sibling pair, is found the samples will appear in `out.fam` without parents. If the direction of the relationship
can't be determined e.g. for parent/child duos a random sample is assigned to be the parent in `out.fam`. The final column
in the `.fam` file specifies how many potential parents the sample had.

##mendel
```

About:   akt mendel - profiles duo/trios
Usage:   ./akt mendel input.bcf -p pedigree.fam

Options:
    -p, --pedigree              pedigree information in plink .fam format
    -o, --out                   output a site only vcf.gz annotated with site specific error rates
    -i, --include               variant filters to apply eg. -i 'TYPE==snp && QUAL>=10 && DP<100000 && HWE<10' 
    -t, --targets [^]<region>   Set regions. Exclude regions with "^" prefix
    -r, --regions <region>      restrict to comma-separated list of regions
    -R, --regions-file <file>   restrict to regions listed in a file
```
* -o : Site only vcf annotated with mendel error rates
* -i : Apply filters see [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions)
* -t : Regions/sites to include/exclude see [bcftools](https://samtools.github.io/bcftools/bcftools.html#common_options)
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -p : The `-p` option specifies the pedigree structure and uses the well-known six column plink .fam format:

1. Family ID
2. Individual ID
3. Paternal ID
4. Maternal ID
5. Sex (1=male; 2=female; other=unknown)
6. Phenotype
Here is an example pedigree file for the [CEPH1463](http://www.nature.com/nmeth/journal/v12/n10/fig_tab/nmeth.3505_SF4.html) pedigree:

```
PG NA12889 0 0  1 -9
PG NA12890 0 0  2 -9
PG NA12891 0 0  1 -9
PG NA12892 0 0  2 -9
PG NA12877 NA12889 NA12890      1 -9
PG NA12878 NA12891 NA12892      2 -9
PG NA12879 NA12877 NA12878      2 -9
PG NA12880 NA12877 NA12878      2 -9
PG NA12881 NA12877 NA12878      2 -9
PG NA12882 NA12877 NA12878      1 -9
PG NA12883 NA12877 NA12878      1 -9
PG NA12884 NA12877 NA12878      1 -9
PG NA12885 NA12877 NA12878      2 -9
PG NA12886 NA12877 NA12878      1 -9
PG NA12887 NA12877 NA12878      2 -9
PG NA12888 NA12877 NA12878      1 -9
PG NA12893 NA12877 NA12878      1 -9
```
##pca
```
Performs principal component analysis on a vcf/bcf
Usage:
./akt pca input.bcf
	 -T --targets-file:		intersecting VCF
	 -o --output:			output vcf
	 -O --outputfmt:		output vcf format
	 -r --regions:			chromosome region
	 -R --regions-file:		list of regions, file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
	 -h --thin:			keep every t variants
	 -m --maf:			minimum MAF
	 -w --weight:			VCF with weights for PCA
	 -N --npca:			first N principle components
	 -a --alg:			exact SVD (slow)
	 -C --covdef:			definition of SVD matrix: 0=(G-mu) 1=(G-mu)/sqrt(p(1-p)) 2=diag-G(2-G) default(1)
	 -e --extra:			extra vectors for Red SVD

```
* -T : Indexed VCF file containing intersecting sites and relevant site info. 
* -o : Output to vcf.
* -O : Output format of vcf b=compressed bcf, z=compressed vcf.
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -s : Comma-separated list of samples to include or exclude if prefixed with "^".
* -S : File of sample names to include or exclude if prefixed with "^". One sample per line.
* -h : Hop this many markers forward from last accepted marker. 
* -m : Minimum MAF a site to be counted. When reading: if AF > 0.5 checks 1-AF > m. 
* -w : Use precalculated principle components.
* -N : Number of principle components to calculate.
* -a : Use JacobiSVD PCA algorithm, which is exact to float precision but very slow.
* -e : Default PCA calculation is the inexact `RedSVD` algorithm, which requires this parameter. The higher the number the more accurate principle components will be obtained. 
* -C : Which matrix to take the PCA of. 0 uses mean subtracted genotype matrix; 1 uses mean subtracted and normalized genotype
matrix; 2 uses normalized covariance matrix with bias term subtracted from diagonal elements.  

`pca` can either calculate principle components from data and do projections onto these components or use precalculated 
principle component vectors
```
~/ancestry_tools$ ./akt pca multisample.bcf -O b -o pca.bcf -m 0.05 -h 10 -N 5 > projections 
```
The file `projections` contains
```
SAMPLE_ID0 P0 P1 P2 P3 P4
SAMPLE_ID1 P0 P1 P2 P3 P4
...
```
The bcf file `pca.bcf` contains
```
bcftools query -f "%INFO/WEIGHT\n" pca.bcf
pc00 pc01 pc02 pc03 pc04
pc10 pc11 pc12 pc13 pc14
...
```
First index is the site index and second is which principle component. If you already have precalculated PCA weights and 
you want to project new samples onto those weights
```
./akt pca new_multisample.bcf -w pca.bcf > projections
```
##cluster
```
Clustering on text files
Usage:   ./akt cluster input.txt
	 -k --K:			number of clusters
	 -i --seed:			random seed for starting values
	 -a --alg:			clustering algorithm 0 = k++means, 1 = gaussian mixture, 2 = density method
	 -c --cols:			column range to read
	 -C --cfile:			initial guess for cluster centres in text file
	 -o --outputcfile:		assigned cluster centres
	 -c --cols:			which columns to use e.g. 2-4
	 -I --maxits:			max number of iterations to use: a=0,1
	 -d --dc:			radius for density method: a=2
	 -p --rho_min:			min density for cluster centre: a=2
	 -D --delta_min:		min radius for cluster centre: a=2
	 --density-plot:		plot the density and finish: a=2
	 -e --silhouette:		calculate silhouette score
```
* -k : Number of clusters. Examine the data to guess this or analyse silhouette scores.
* -i : Random seed.
* -a : 0 to use k++ means clustering, 1 to use EM with Gaussians, 2 to use density method.
* -C : file with initial guess for cluster centres. If analysing 2d data this should contain 2 columns and K rows.
* -c : Which columns in input file to use.
* -I : Maximum number of iterations to use for alg 0 or 1.
* -d : radius around each point for counting density.
* -p : min density for cluster centre.
* -D : min radius for cluster centre.
* --density-plot : plot delta-density graph and finish.
* -e : Calculate silhouette score (goodness of cluster assignment).

e.g. cluster the first three principle components from the output of vcfpca using k++ means. 
Different clusters should correspond to different ancestries
```
./akt cluster projections -k 4 -c 2-4 -e > clustered
```
Each block in `clustered` contains all the data for points in that cluster with the row layout
```
P2 P3 P4 (SILHOUETTE) CLUSTERID P1 P5 ...
```
To visually inspect the clustered data try the following gnuplot command
```
gnuplot> splot for [i=0:3:1] 'clustered' every:::i::i
```
The silhouette column can be used to filter on badly clustered samples. 
~1 means a sample is closer to the centroid of its assigned cluster than the next nearest cluster.
A sample can be badly clustered for many reasons: bad choice of K, convergence of the clustering algorithm 
to bad local minimum, individual is of mixed ancestry. You could filter out badly clustered samples using
```
awk '{ if($4 > 0.5) print $0}' clustered > well_clustered
```
###density clustering
This is a better clustering algorithm but requires some tuning by hand. First run
```
./akt cluster projections -c 2-4 -a 2 -d 1 --density-plot > dplot
```
the file dplot is a two column data file of density versus delta 
(the distance to the closest data point of higher density). By examining this plot outliers 
are obvious and a cutoff in density and delta can be determined that isolates these points. 
Once this cutoff is determined run
```
./akt cluster projections -c 2-4 -a 2 -d 1 -p 0 -D 75 > clustered
```
Where 0 and 75 are the cutoffs determined from `dplot`. This will output the number of clusters that were found 
(including one for unassigned points) and `clustered` will contain the clustered data in the usual format.
```
gnuplot> splot for [i=0:4:1] 'clustered' every:::i::i
```
Unassigned data points are always put in cluster 0. The silhouette score for unassigned data is undefined.

##stats
```
Calculate correlation across a population
Usage:
./akt stats input_filename.vcf
Expects input_filename.vcf to contain hard genotypes
	 -F --flank:			size of left and right flanking region (1000bp)
	 -b --block:			-F argument is number of markers instead of number of base pairs
	 -x --afonly:			calculate allele freq only
	 -c --output_cor:		output sitewise correlation (false)
	 -C --output_cormin:		output sitewise correlation greater than (0)
	 -a --aftag:			allele frequency tag
	 -o --output:			output vcf
	 -O --outputfmt:		output vcf format
	 -r --regions:			chromosome region
	 -R --regions-file:		list of regions, file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
```
* -F : Correlation with variants in a window of size f base pairs to left and right or each variant. 
* -b : The number in the -F argument now interpreted as number of flanking variants instead of flanking positions.
* -x : Do not calculate covariance matrix.
* -c : Output all correlation.
* -C : Output correlation values greater than this.
* -a : Allele frequency tag e.g. "TAG" reads or writes "TAG_AF" in the target VCF. 
* -o : Output to vcf.
* -O : Output format of vcf b=compressed bcf, z=compressed vcf.
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -s : Comma-separated list of samples to include or exclude if prefixed with "^".
* -S : File of sample names to include or exclude if prefixed with "^". One sample per line.

This tool lets us calculate allele frequencies and correlation matrices from multisample vcfs. Running 
```
./akt stats panel.bcf -F 10000 -O b -o tmp.bcf -a "TAG" -c
```
The output, tmp.bcf, has no sample columns and info fields 
* TAG_AF : Allele frequencies.
* TAG_SIG : Allele variances.
* TAG_LD : LDscore = sum of COR.
* TAG_COR : correlation of this site with neighbouring sites.
* TAG_CORP : position of neighbouring sites.

Correlation is bounded by [-1,1]. With the `-f 10000` option it will output the correlation of each site with the 
sites closer than 10000bp to the left and right. The LD metric is described in [this paper](http://www.nature.com/ng/journal/v47/n3/full/ng.3211.html)

##ibd
```
Find IBD regions from a VCF
Usage:
./akt ibd in.bcf -p sites.bcf
Expects input_filename.vcf to contain hard genotypes
	 -T --targets-file:		intersecting VCF
	 -r --regions:			chromosome region
	 -R --regions-file:		list of regions, file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
	 -h --thin:			keep every t variants
	 -a --aftag:			allele frequency tag
	 -f --pairfile:			file containing sample pairs
	 -n --nthreads: 		num threads
	 -m --maf:			minimum MAF
	 -e --error: 			If no GQ then this is error probability(1e-3)
	 -M --thresh: 			likelihood output threshold (default 50)
	 -L --length: 			length output threshold (default 100000)
	 -w --wsize: 			window size (default 20)
	 -x --maxerr: 			stop when encountering a window with this many ( 0/0 , 1/1 ) pairs (default 2)
	 -l --lsize: 			try to join long regions closer than this (default 100000)
```
* -T : Indexed VCF file containing intersecting sites and relevant site info. 
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -s : Comma-separated list of samples to include or exclude if prefixed with "^".
* -S : File of sample names to include or exclude if prefixed with "^". One sample per line.
* -h : Hop this many markers forward from last accepted marker. 
* -a : Allele frequency tag e.g. "TAG" reads or writes "TAG_AF" in the target VCF. 
* -f : File containing pairs of samples to compute, 2 column: Sample1 Sample2. 
* -n : Number of threads to use. 
* -m : Minimum MAF a site to be counted. When reading: if AF > 0.5 checks 1-AF > m. 
* -e : Probability of error per base
* -M : Lower bound of likelihood ratio score
* -L : Lower bound of segment length
* -w : Number of markers per window
* -x : Number of hom-ref hom-alt matches allowed in a window before breaking a region in first stage.
* -l : Minimum IBD region to output

This tool lets us find regions which are shared IBD between two samples. 
```
./akt ibd input.bcf -p panel.bcf -r 2 > sharing
```
The file `sharing` is in the format
```
SAMPLE1 SAMPLE2 START END SCORE MATCHES LEN
```
Where higher scores indicate a segment that is more likely to be shared IBD between samples 1 and 2. Matches
are the number of IBS consistent matches between the two segments and LEN=END-START is the length of the segment.

The algorithm is a variation of [EAGLE](https://data.broadinstitute.org/alkesgroup/Eagle/). First we scan along the
genome in windows of `w` counting the number of opposite homozygotes. With no errors an opposite homozygote
indicates that the two windows cannot be IBD. Once we have counted in all the windows we join up neighbouring 
windows to create seed regions and attempt to extend these regions through windows with (w-x-1) opposite homozygotes. 
We calculate a running score based on the method in [this paper](http://biorxiv.org/content/early/2015/12/18/028282). 
Once these segments are calculated we try to join nearby segments (closer than `l`) if the score of the combined
segment is greater than the score of either one separately. By filtering on segment length and score this code does 
a reasonable job of identifying IBD shared segments longer than about 1Mb from unphased data.

##LDplot
```
Calculate correlation across a population
Usage:
./akt ldplot input_filename.vcf
Expects input_filename.vcf to contain hard genotypes
	 -r --regions:			chromosome region
	 -R --regions-file:		list of regions, file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
```
* -r : Comma-separated list of regions, chr:from-to.
* -R : File containing 3 columns: CHROM, POS and POS_TO. 
* -s : Comma-separated list of samples to include or exclude if prefixed with "^".
* -S : File of sample names to include or exclude if prefixed with "^". One sample per line.
This generates the correlation matrix of all included variants.
```
./akt LDplot input.bcf -r 20:0-200000 > sigma
```
the file `sigma` contains the correlation matrix of all variants in the region `20:0-200000` and can be plotted
with various tools e.g.
```
gnuplot> plot 'sigma' matrix with image
```
##admix
```
Approximate admixture fractions
Usage:   ./akt admix input.txt -C centres -c 2-3 
	 -c --cols:			column range to read
	 -C --cfile:			initial guess for cluster centres in text file

```
* -c : Which columns in input file to use.
* -C : file with vectors whose transformation to admixture fractions is known.

`admix` attempts to use the output of `pca` to assign admixture fractions to data based on known populations.
The input file should be in the form
```
e1 e2 e3 ... ed POP a1 a2 a3 ... ad
...
```
where `e1...` are coordinates in the space of principal components and `a1...` are coordinates in the space
of admixture fractions. The mapping from PCA space to admixture space is
```
Ta + v = e
```
with inverse
```
a = T^{-1}(e-v)
```
`admix` calculates this inverse mapping.

```
~/ancestry_tools$ ./akt admix pcadata -c 2-6 -C pca_to_admix.txt > admixtures
```
The file `admixtures` contains
```
SAMPLE_ID %POP1 %POP2 ... %POPD
...
```
A vector must be chosen to fix the zero point of the transformation (the vector `v` above). 
For K populations the `-C` input file can contain K or K+1 lines.
If K+1 vectors are given the one with the lowest norm in the space of admixtures equals `v`.
If K vectors and their images in admixture space are specified then the zero point is fixed by 
choosing an arbitrary vector orthogonal to the sides of the tetrahedron made by the K inputs
and roughly equidistant from them. There are two possible choices of sign: the positive one is always chosen.



##Example Workflow
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
./akt kin ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R data/1000G.snps.nochr.vcf.gz -n 4 > test_kin
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
##Segment sharing
We can try to find long regions shared IBD between two samples. 
We choose to demonstrate this for samples which are closely related, though not flagged by the kin tool or by 1000 Genomes. 
These are all East Asian individuals. 
```
cat pairs.txt
HG00464 HG00537
HG00477 HG00544
HG00477 HG00542
HG00475 HG00544
HG00475 HG00542
```
Then run
```
./akt ibd ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -T sites.EAS.bcf -a "EAS" -f pairs.txt > test_ibd
```
We only look at chromosome 20 and a small number of pairs. `akt ibd` uses the LD score metric computed by `akt stats`. If
this is not present the markers are assumed to be in linkage equilibrium. Unfortunately there isn't good truth data for IBD. 
We run  [IBDseq](http://faculty.washington.edu/browning/ibdseq.html) on the same data
```
bcftools view ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -S test.Cluster1 -O z -o test.20.vcf.gz;
java -jar ibdseq.r1206.jar gt=EAS.test.20.vcf.gz out=Seq.EAS.test.20
```
Comparing the IBD segments gives

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/segment.png)

Where the yscale is arbitrary and the xscale is position along chromosome 20. We show the segments 
as reported by IBDseq in red and by akt in black. For each of the 5 pairs examined 
(chosen because of the presence of significant IBD sharing) the two methods give similar answers. 
We expect IBDseq to be more accurate since it uses phased haplotypes where akt uses only genotypes.
Compared to IBDseq, akt produces some false positives and some broken segments reported as continuous 
by IBDSeq. These can be controlled by tuning the score and minimum segment length thresholds.
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
`-w` option and the file `1000G.snps.nochr.vcf.gz` or `1000G.snps.withchr.vcf.gz`
```
./akt admix test_pcaproj -c 2-6 -C data/1000G.pca_to_admix > test_alladmix
```

![alt text](https://raw.githubusercontent.com/Illumina/akt/master/docs/test_alladmix.png)


