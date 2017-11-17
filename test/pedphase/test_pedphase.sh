#!/usr/bin/env bash

for i in pedphase/test*.vcf.gz;
do
    echo Testing $i
    ../akt pedphase $i -o - 2> /dev/null  | grep -A1000 CHROM > ${i%vcf.gz}observed
    diff ${i%vcf.gz}observed ${i%vcf.gz}expected
    ../akt pedphase $i -o - 2> /dev/null | grep -A1000 CHROM  > ${i%vcf.gz}observed
    diff ${i%vcf.gz}observed ${i%vcf.gz}expected
    echo "PASSED"
done

## how to update tests - USE WITH CAUTION
# for i in pedphase/test*.vcf.gz;
# do
#     echo Testing $i
#     ../akt pedphase $i -p pedphase/pedigree.fam -o - 2> /dev/null  | grep -A1000 CHROM > ${i%vcf.gz}expected
# done
