#Using akt

akt uses the syntax
```
./akt COMMAND <options>
```
To see a list of available tools use
```
replace_akt_run
```

##common options 
There are a number of options that are shared by multiple akt subcommands which we list here. We have tried to keep these consistent with [bcftools](http://samtools.github.io/bcftools/bcftools.html) where possible.

**-R** *FILE* a file (tabixed VCF or bed) containing the markers to perform analysis on. **-R**/**-r** uses tabixes jumping for fast look up  
**-r** *REGION* same as **-R** but a string containing the region eg. `chr1:1000000-2000000`  
**-T** *FILE* same as **-R** but streams rather than tabix jumps ie. is slow  
**-t** *TARGET* same as **-r** but streams rather than tabix jumps ie. is slow  
**-S** *SAMPLES* File of sample names to include or exclude if prefixed with "^"  
**-s** *SAMPLES* Comma-separated list of samples to include or exclude if prefixed with "^"  
**-n** *VALUE* Number of threads to use.  
**-o** *FILE* Output file name  
**-O** *v/z/b/u* Output format of vcf b=compressed bcf, z=compressed vcf, u=uncompressed bcf, v=uncompressed vcf
**-m** *VALUE* Minimum MAF a site to be counted

