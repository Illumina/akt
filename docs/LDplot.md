##LDplot

This generates the correlation matrix of all included variants.

regions_option
regionsfile_option
samples_option
samplesfile_option

```
./akt LDplot input.bcf -r 20:0-200000 > sigma
```

the file `sigma` contains the correlation matrix of all variants in the region `20:0-200000` and can be plotted
with various tools e.g.

```
gnuplot> plot 'sigma' matrix with image
```
