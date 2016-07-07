##thin

Performs tag SNP selection . This is useful for choosing a subset of variants to use in analyses that do not require anywhere near the number of markers typically obtained from WGS (notably `pca`/`kin`) and assume linkage equilibrium between markers. Note this routine is rather different to PLINK's pruning method. We attempt to find the *K* (non-redundant) variants that tag as much variation as possible via a greedy algorithm. Whereas PLINK prunes away variants that can be predicted from other variants. In practice, for applications such as PCA/kinship calculations, it appears any set of reasonably common and sparse markers is appropriate.

The default values are appropriate for WGS data after filtering for MAF>=5% (about 7 million SNPs).

output_option
**-w** *VALUE* window size within to consider possible tag SNPs (default 250)
**-c** *VALUE* minimum r^2 value to consider a SNP tagged (default 0.5)
**-k** *VALUE* number of tag SNPs to selection (default 20000)
nthread_option

###example
Here was filter low frequency variants (and indels) with bcftools and pipe the output straight to `akt thin` to select a subset of markers:

```
 bcftools view -t ^X -v snps -i 'MAF>=0.05 && N_ALT==1' ~/kimura/resources/1000G/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -Ou  | ./akt thin -n 48 - | bgzip -c > thinned.vcf.gz
```

**Note:** This routine is somewhat slow (the above command takes a few hours). We provided a set of reliable thinned SNPs for both WGS and Exome data on builds 37 and 38 of the human genome. These site-only VCFs are under `data/` in the akt repository.

