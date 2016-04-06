##range
```
replace_range_run
```
* regions_option
* regionsfile_option
This tool takes an input VCF and outputs a regions-file.
```
./akt range input.bcf > regions
```
the file `regions` contains
```
CHROM1 POS POS_TO
CHROM2 POS POS_TO
...
```
where `CHROM1` etc. are the distinct chromosomes found in `input.bcf`, `POS` is the position of the first
variant of that chromosome and `POS_TO` is the position of the last variant on that chromosome. Among other
things this is useful for running the pca or kin tools with the `--enable-pi` option.
