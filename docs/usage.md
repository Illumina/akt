#Using akt

akt uses the syntax
```
./akt COMMAND <options>
```
To see a list of available tools use
```

Program:	akt (Ancestry and Kinship Tools)
Version:	45c10fa
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
	metafreq                 examine two files for AF differences

```

##common options 
There are a number of options that are shared by multiple akt subcommands which we list here. We have tried to keep these consistent with [bcftools](http://samtools.github.io/bcftools/bcftools.html) where possible.

**-R** *FILE* a file (tabixed VCF or bed) containing the markers to perform analysis on. **-R**/**-r** uses tabixes jumping for fast look up  
**-r** *REGION* same as **-R** but a string containing the region eg. `chr1:1000000-2000000`  
**-T** *FILE* same as **-R** but streams rather than tabix jumps ie. is slow  
**-t** *TARGET* same as **-r** but streams rather than tabix jumps ie. is slow  
**-S** *SAMPLES* File of sample names to include or exclude if prefixed with "^"  
**-s** *SAMPLES* Comma-separated list of samples to include or exclude if prefixed with "^"  
**-n** *VALUE* Number of threads to use.  
**-o** *FILE* Output file name
**-O** *v/z/b/u* Output format of vcf b=compressed bcf, z=compressed vcf, u=uncompressed bcf, v=uncompressed vcf
**-m** *VALUE* Minimum MAF a site to be counted

##pca

Performs principal component analysis on a BCF/VCF. Can also be used to project samples onto pre-calculated principal components from another cohort. Uses a randomised SVD by default for very fast computation. WGS data is far denser than required for a meaningful PCA, it is recommended you provide a thinned set of sites via the `-R` command.

**-o** *FILE* see common options  
**-O** *z|b|v|u* see common options  
**-r** *REGION* see common options  
**-R** *FILE* see common options  
**-T** *FILE* see common options  
**-s** *SAMPLES* see common options  
**-S** *SAMPLES* see common options  
**-W** *FILE* Use precalculated principle components.
**-N** *VALUE* Number of principle components to calculate.
**-a** Use JacobiSVD PCA algorithm, which is exact to float precision but very slow.
**-e** Default PCA calculation is the inexact `RedSVD` algorithm, which requires this parameter. The higher the number the more accurate principle components will be obtained. 
**-F** File to output the singular values.
**-C** Which matrix to take the PCA of. 0 uses mean subtracted genotype matrix; 1 uses mean subtracted and normalized genotype
matrix; 2 uses normalized covariance matrix with bias term subtracted from diagonal elements.  


```
./akt pca multisample.bcf -R data/wgs.grch37.vcf.gz -O b -o pca.bcf > pca.txt
```

The file `pca.txt` contains
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
First index is the site index and second which is the coefficient (loading) that can be used to project other samples onto these principal components. For example we could project a new set of samples onto these same PCs via:
```
./akt pca new_multisample.bcf -W pca.bcf > projections
```
##kin

**-k** *value* Only output pairs with kinship coefficient greater than *value*  
**-F** *FILE* a file containing population allele frequencies to use in kinship calculation  
**-M** *0/1/2* type of estimator.  0:[plink (default)](https://www.cog-genomics.org/plink2/ibd) 1:[king-robust](http://bioinformatics.oxfordjournals.org/content/26/22/2867.full) 2:[genetic-relationship-matrix](http://cnsgenomics.com/software/gcta/estimate_grm.html)  
**-a** *TAG* Allele frequency tag e.g. "TAG" reads or writes "TAG_AF" in the target VCF.  
**-n** *VALUE see common options  
**-R** *FILE* see common options  
**-r** *REGION* see common options  
**-T** *FILE* see common options  
**-t** *TARGET* see common options  
**-S** *SAMPLES* see common options  
**-s** *SAMPLES* see common options  


Run the kinship calculation by giving akt a multi-sample vcf/bcf file:

```
$ akt kin multisample.bcf -R data/wgs.grch37.vcf.gz -n 32 > kin.txt
```
###Choice of estimator


The default algorithm (`-M 0`) used to calculate IBD is taken from [PLINK](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/) with some minor changes. It outputs the following seven column format:

```
ID1 ID2 IBD0 IBD1 IBD2 KINSHIP NSNP
```

As with PLINK, we set KINSHIP = 0.5 * IBD2 + 0.25 * IBD1. Our IBD values may slighly differ to PLINK's (by desing) due to the following differences:

* No 'bias correction' since allele frequencies are assumed to be accurate
* Normalization as follows:
	* if IBD0 > 1: IBD0 = 1, IBD1 = 0, IBD2 = 0
	* if IBD1 < 0: IBD1 = 0
	* if IBD2 < 0: IBD2 = 0
	* norm = IBD0 + IBD1 + IBD2
	* IBD0 /= norm, IBD1 /= norm; IBD2 /= norm;
* We do **not** follow PLINK which forces IBD to obey consistency conditions - this affects the clustering that is required for the `relatives` code.

The second method (`-M 1`) uses the robust kinship coefficent estimate describing in the [KING paper](http://bioinformatics.oxfordjournals.org/content/26/22/2867.full). This may be preferable when your cohort has large amounts of population structure. The IBD estimates and output format are the same as for `-M 0`.

The third method (`-M 2`) outputs the genetic relationship matrix used by [GCTA](http://cnsgenomics.com/software/gcta/) (among other software).  The output format is `ID1 ID2 GR` where `GR` is the genetic correlation between samples `ID1` and `ID2`. This matrix is important for fitting linear mixed-effect models but we have found methods 0 and 1 more appropriate for sample QC purposes.

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
* **-r** *REGION* see common options  
* **-R** *FILE* see common options  
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
	 -a --aftag:			allele frequency tag (default AF)
	 -o --output:			output vcf
	 -O --outputfmt:		output vcf format
	 -r --regions:			chromosome region
	 -R --regions-file:		restrict to regions listed in a file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
```
* -F : Correlation with variants in a window of size f base pairs to left and right or each variant. 
* -b : The number in the -F argument now interpreted as number of flanking variants instead of flanking positions.
* -x : Do not calculate covariance matrix.
* -c : Output all correlation.
* -C : Output correlation values greater than this.
* aftag_option
* **-o** *FILE* see common options  
* **-O** *z|b|v|u* see common options  
* **-r** *REGION* see common options  
* **-R** *FILE* see common options  
* **-s** *SAMPLES* see common options  
* **-S** *SAMPLES* see common options  

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
	 -T --targets-file:		similar to -R but streams rather than index-jumps
	 -r --regions:			chromosome region
	 -R --regions-file:		restrict to regions listed in a file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
	 -h --thin:			keep every h variants
	 -a --aftag:			allele frequency tag (default AF)
	 -f --pairfile:			file containing sample pairs to perform calculations on
	 -n --nthreads: 		num threads
	 -m --maf:			minimum MAF
	 -e --error: 			If no GQ then this is error probability(1e-3)
	 -M --thresh: 			likelihood output threshold (default 50)
	 -L --length: 			length output threshold (default 100000)
	 -w --wsize: 			window size (default 20)
	 -x --maxerr: 			stop when encountering a window with this many ( 0/0 , 1/1 ) pairs (default 2)
	 -l --lsize: 			try to join long regions closer than this (default 100000)
```
* **-T** *FILE* see common options  
* **-r** *REGION* see common options  
* **-s** *SAMPLES* see common options  
* **-S** *SAMPLES* see common options  
* thin_option
* aftag_option
* pairfile_option
* **-n** *VALUE* see common options  
* **-m** *VALUE* see common options  
* -e : Probability of error per base
* -M : Lower bound of likelihood ratio score
* -L : Lower bound of segment length
* -w : Number of markers per window
* -x : Number of hom-ref hom-alt matches allowed in a window before breaking a region in first stage.
* -l : Minimum IBD region to output

This tool lets us find regions which are shared IBD between two samples. 
```
./akt ibd input.bcf -R panel.bcf -r 2 > sharing
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
	 -R --regions-file:		restrict to regions listed in a file
	 -s --samples:			list of samples
	 -S --samples-file:		list of samples, file
```
* **-r** *REGION* see common options  
* **-R** *FILE* see common options  
* **-s** *SAMPLES* see common options  
* **-S** *SAMPLES* see common options  
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



##metafreq
```
Compare AFs between two cohorts
Usage:
akt metafreq a.vcf.gz b.vcf.gz -Oz -o meta.sites.vcf.gz
	 -r --regions:			chromosome region
	 -R --regions-file:		restrict to regions listed in a file
	 -T --targets-file:		similar to -R but streams rather than index-jumps
	 -o --output:			output vcf
	 -O --outputfmt:		output vcf format
	 -a --aftag:			allele frequency tag (default AF)
```
* **-T** *FILE* see common options  
* **-r** *REGION* see common options  
* **-R** *FILE* see common options  
* aftag_option
* **-o** *FILE* see common options  
* **-O** *z|b|v|u* see common options  


This tool lets us find sites with statistically significant differences in frequency between two samples. 
```
./akt ibd f1.bcf -p f2.bcf -O b -o out.bcf 
```
The file `out.bcf` contains two fields `QF` and `QX` equal to -log10 of the p-value
calculated using Fisher's exact test and the chi-squared test, respectively.
Higher scores indicate a more statistically significant difference. These values are capped at 100
and typically the chi-squared test is more likely to reject the null-hypothesis that
the two cohorts have identical distributions.

