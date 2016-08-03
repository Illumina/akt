##prune

This is similar to PLINK's pairwise LD-pruning routine. The algorithm slides along the genome, and calculates the squared correlation coefficient (r2) between plus and minus *b* flanking variants of a given variant. If r2 is greater than the specified threshold, than the variant with the **higher** MAF is removed. For WGS data, this can still result in a very dense set of markers due to the large number of rare variants that are in LD with very few other markers.


Example:
```
 bcftools view -v snps -i 'MAC>1' -Ou uk10.chr20.bcf | akt prune - -Ou | bcftools view -Ob -o uk10k.chr20.pruned.bcf
```
