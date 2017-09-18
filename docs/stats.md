## stats

This tool lets us calculate allele frequencies and correlation matrices from multisample vcfs.

**-F** *VALUE* Correlation with variants in a window of size f base pairs to left and right or each variant.  
**-b** The number in the -F argument now interpreted as number of flanking variants instead of flanking positions.  
**-x** Do not calculate covariance matrix.  
**-c** Output all correlation.  
**-C** Output correlation values greater than this.  
aftag_option
output_option
outputfmt_option
regions_option
regionsfile_option
samples_option
samplesfile_option

 Running 
```
./akt stats panel.bcf -F 10000 -O b -o tmp.bcf -a "TAG" -c
```
The output, tmp.bcf, has no sample columns and info fields 
* TAG_AF : Allele frequencies.
* TAG_SIG : Allele variances.
* TAG_LD : LDscore = sum of COR.
* TAG_COR : correlation of this site with neighbouring sites.
* TAG_CORP : position of neighbouring sites.

Correlation is bounded by [-1,1]. With the `-f 10000` option it will output the correlation of each site with the 
sites closer than 10000bp to the left and right. The LD metric is described in [this paper](http://www.nature.com/ng/journal/v47/n3/full/ng.3211.html)

