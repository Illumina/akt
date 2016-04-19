##metafreq
```
replace_metafreq_run
```
* targetfile_option
* regions_option
* regionsfile_option
* aftag_option
* output_option
* outputfmt_option


This tool lets us find sites with statistically significant differences in frequency between two samples. 
```
./akt ibd f1.bcf -p f2.bcf -O b -o out.bcf 
```
The file `out.bcf` contains two fields `QF` and `QX` equal to -log10 of the p-value
calculated using Fisher's exact test and the chi-squared test, respectively.
Higher scores indicate a more statistically significant difference. These values are capped at 100
and typically the chi-squared test is more likely to reject the null-hypothesis that
the two cohorts have identical distributions.

