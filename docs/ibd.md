##ibd

Detects large chunks of DNA shared Identical-by-Descent (IBD) between pairs of samples.

targetfile_option
regions_option
samples_option
samplesfile_option
nthread_option
maf_option
**-h** *VALUE* keep every *VALUE*th marker  
**-a** *STRING* allele frequency tag (default AF)  
**-f** *FILE* file containing sample pairs to perform calculations on  
**-e** *VALUE* Probability of error per base  
**-M** *VALUE* Lower bound of likelihood ratio score  
**-L** *VALUE* Lower bound of segment length  
**-w** *VALUE* Number of markers per window  
**-x** *VALUE* Number of hom-ref hom-alt matches allowed in a window before breaking a region in first stage.  
**-l** *VALUE* Minimum IBD region to output  

This tool lets us find regions which are shared IBD between two samples. 
```
./akt ibd input.bcf -R panel.bcf -r 2 > sharing
```
The file `sharing` is in the format
```
SAMPLE1 SAMPLE2 START END SCORE MATCHES LEN
```
Where higher scores indicate a segment that is more likely to be shared IBD between samples 1 and 2. Matches
are the number of IBS consistent matches between the two segments and LEN=END-START is the length of the segment.

The algorithm is a variation of [EAGLE](https://data.broadinstitute.org/alkesgroup/Eagle/). First we scan along the
genome in windows of `w` counting the number of opposite homozygotes. With no errors an opposite homozygote
indicates that the two windows cannot be IBD. Once we have counted in all the windows we join up neighbouring 
windows to create seed regions and attempt to extend these regions through windows with (w-x-1) opposite homozygotes. 
We calculate a running score based on the method in [this paper](http://biorxiv.org/content/early/2015/12/18/028282). 
Once these segments are calculated we try to join nearby segments (closer than `l`) if the score of the combined
segment is greater than the score of either one separately. By filtering on segment length and score this code does 
a reasonable job of identifying IBD shared segments longer than about 1Mb from unphased data.

