#!/usr/bin/env bash

./akt pedphase test/pedphase/test1.vcf.gz -p test/pedphase/pedigree.fam -o - | bcftools view -H > test/pedphase/test1.out
diff test/pedphase/test1.out test/pedphase/test1.vcf

./akt pedphase test/pedphase/test2.vcf.gz -p test/pedphase/pedigree.fam -o - | bcftools view -H > test/pedphase/test2.out
diff test/pedphase/test2.out test/pedphase/test2.vcf

