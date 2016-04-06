##pca
```
replace_pca_run

```
* targetfile_option
* output_option
* outputfmt_option
* regions_option
* regionsfile_option
* samples_option
* samplesfile_option
* thin_option
* maf_option
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
