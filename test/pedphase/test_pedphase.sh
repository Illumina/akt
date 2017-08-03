#!/usr/bin/env bash

echo "Test 1"
./akt pedphase test/pedphase/test1.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test1.out
diff test/pedphase/test1.out test/pedphase/test1.vcf
echo "PASSED"

echo "Test 2"
./akt pedphase test/pedphase/test2.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test2.out
diff test/pedphase/test2.out test/pedphase/test2.vcf
echo "PASSED"

echo "ALL TESTS PASSED!"
