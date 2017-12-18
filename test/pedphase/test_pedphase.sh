#!/usr/bin/env bash

##small smoke tests
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
#      ../akt pedphase $i -o - 2> /dev/null  | grep -A1000 CHROM > ${i%vcf.gz}expected
# done


##a chr20 test for accuracy
# ../akt pedphase pedphase/Pedigree.chr20.vcf.gz | bcftools view -i 'PS=="."' -Oz -o /tmp/phased.vcf.gz
# vcftools --gzvcf /tmp/phased.vcf.gz --gzdiff pedphase/NA12878.truth.grch37.chr20.vcf.gz --diff-switch-error
# grep NA12878 out.diff.indv.switch
