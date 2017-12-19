# akt - ancestry and kinship toolkit

Copyright (c) 2017, Illumina, Inc. All rights reserved. This software is not commercially supported.

Ancestry and Kinship Tools (AKT) provides a handful of useful statistical genetics routines using the [htslib](http://www.htslib.org/) API for input/output. This means it can seamlessly read BCF/VCF files and play nicely with [bcftools](http://samtools.github.io/bcftools/bcftools.html).

Please cite the [AKT pre-print](http://biorxiv.org/content/early/2017/04/10/047829) if you find this software useful.

## License and dependencies

AKT is freely available under the [GPL3 license](https://github.com/Illumina/agg/blob/master/LICENSE). 

AKT relies on HTSlib and Eigen. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a header-only library for matrix algebra released under the [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/). [HTSlib](http://www.htslib.org/) is a library for efficently parsing vcf/bcf files released under the [MIT/Expat License](http://choosealicense.com/licenses/mit/). Both libraries are included with AKT.

## Installation instructions

```
git clone https://github.com/Illumina/akt.git
cd akt/
make
```
If you get a warning about 'omp.h not found' (especially on OSX)
you can try
```
make no_omp
```
Everything will be run on a single thread, so the `-n` option does nothing
in `akt kin` and `akt ibd`.

## Quick start
akt uses the syntax
```
./akt COMMAND <options>
```

To see a list of available functionality use:
```
Program:        akt (Ancestry and Kinship Tools)
Version:        v0.3.9-33-g2f4bfea
Copyright (c) 2017, Illumina, Inc. All rights reserved. See LICENSE for further details.

Usage:  akt <command> [options]

        pca                      principal component analysis
        kin                      calculate kinship coefficients
        relatives                discover pedigrees
        unrelated                generate a list of unrelated individuals
        pedphase                 Mendelian transmission phasing for duos/trios
```

Some useful routine analyses (assumes you are in the akt directory):
```
#Do a PCA:
./akt pca -R data/wgs.grch37.vcf.gz input.bcf > pca.txt
Rscript scripts/pca.R pca.txt

#Project samples onto 1000G principal components
./akt pca -W data/wgs.grch37.vcf.gz input.bcf > 1000G.pca.txt
Rscript scripts/1000G_pca.R 1000G.pca.txt

#Calculate some (robust) kinship coefficients
./akt kin -R data/wgs.grch37.vcf.gz -M 1 input.bcf > kinship.txt

#give me duplicated sample pairs:
awk '{if($6>.4) print $0}' kinship.txt

#reconstruct pedigrees:
akt relatives kinship.txt -p pedigree

#give me a list of unrelated samples
akt unrelated kinship.txt > unrelated.ids

```
Full documentation is available at [http://illumina.github.io/akt](http://illumina.github.io/akt).

 