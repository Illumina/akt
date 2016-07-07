##mendel

Profiles Mendelian inheritance/inconsistency patterns in a sample. Useful for evaluating different sets of filters.

output_option
target_option
regions_option
regionsfile_option
**-i** *EXPRESSION* Apply filters see [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions)
**-p** *FILE* The `-p` option specifies the pedigree structure and uses the well-known six column plink .fam format:

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
