# akt - ancestry and kinship toolkit

Copyright (c) 2016, Illumina, Inc. All rights reserved. This software is not commercially supported.

Ancestry and Kinship Tools (AKT) provides a number of useful statistical genetics routines using the [htslib](http://www.htslib.org/) API for input/output. This means it can seamlessly read BCF/VCF files and play nicely with [bcftools](http://samtools.github.io/bcftools/bcftools.html).

Please cite the [AKT pre-print](http://biorxiv.org/content/early/2016/04/10/047829) if you find this software useful.

##License and dependencies

AKT is freely available under the [GPL3 license](https://github.com/Illumina/agg/blob/master/LICENSE). 

AKT relies on HTSlib and Eigen. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a header-only library for matrix algebra released under the [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/). [HTSlib](http://www.htslib.org/) is a library for efficently parsing vcf/bcf files released under the [MIT/Expat License](http://choosealicense.com/licenses/mit/). Both libraries are included with AKT.

##Installation instructions

```
git clone https://github.com/Illumina/akt.git
cd akt/
make
```
If you get a warning about 'omp.h not found' (especially on Macs)
you can try
```
make no_omp
```
Everything will be run on a single thread, so the `-n` option does nothing
in `akt kin` and `akt ibd`.

##Quick start
akt uses the syntax
```
./akt COMMAND <options>
```

To see a list of available functionality use:
```
Program:	akt (Ancestry and Kinship Tools)
Version:	a963bf1
Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.

Usage:	akt <command> [options]

	pca                      principal component analysis
	kin                      detect average IBD sharing
	relatives                discover pedigrees
	mendel                   profile Mendelian inhertiance and inconsistencies in known pedigrees
	cluster                  perform cluster analyses
	LDplot                   output correlation matrix
	stats                    calculate AF and LD metrics
	prune                     perorms LD pruning of variants
	tag                      selects a set of K tagging variants
	metafreq                 examine two files for AF differences
```

Some useful routine analyses (assumes you are in the akt directory):
```
#Do a PCA:
./akt pca -R data/wgs.grch37.vcf.gz input.bcf > pca.txt
Rscript scripts/pca.R pca.txt

#Project samples onto 1000G principal components
./akt pca -W data/wgs.grch37.vcf.gz input.bcf > 1000G.pca.txt
Rscript scripts/1000G_pca.R 1000G.pca.txt

#Find close relatives/duplicates:
./akt pca -R data/wgs.grch37.vcf.gz input.bcf > kinship.txt
awk '{if($5>.15) print $0}' kinship.txt

#reconstruct pedigrees:
akt relatives -p pedigree kinship.txt 
```
Full documentation is available in [docs/usage.md](docs/usage.md) along with a large set of example use cases in [docs/examples.md](docs/examples.md)

