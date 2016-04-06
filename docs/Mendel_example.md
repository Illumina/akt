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
