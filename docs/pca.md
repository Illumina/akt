##pca

Performs principal component analysis on a BCF/VCF. Can also be used to project samples onto pre-calculated principal components from another cohort. Uses a randomised SVD by default for very fast computation. WGS data is far denser than required for a meaningful PCA, it is recommended you provide a thinned set of sites via the `-R` command.

output_option
outputfmt_option
regions_option
regionsfile_option
targetfile_option
samples_option
samplesfile_option
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
