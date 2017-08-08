#!/usr/bin/env bash

echo "Test 1"
./akt pedphase test/pedphase/test1.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test1.out
diff test/pedphase/test1.out test/pedphase/test1.vcf
./akt pedphase test/pedphase/test1.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test1.out
diff test/pedphase/test1.out test/pedphase/test1.vcf
echo "PASSED"

echo "Test 2"
./akt pedphase test/pedphase/test2.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test2.out
diff test/pedphase/test2.out test/pedphase/test2.vcf
./akt pedphase test/pedphase/test2.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test2.out
diff test/pedphase/test2.out test/pedphase/test2.vcf
echo "PASSED"

echo "Test 3"
./akt pedphase test/pedphase/test3.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test3.out
diff test/pedphase/test3.out test/pedphase/test3.vcf
./akt pedphase test/pedphase/test3.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test3.out
diff test/pedphase/test3.out test/pedphase/test3.vcf
echo "PASSED"

echo "Test 5"
./akt pedphase test/pedphase/test5.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test5.out
diff test/pedphase/test5.out test/pedphase/test5.vcf
./akt pedphase test/pedphase/test5.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test5.out
diff test/pedphase/test5.out test/pedphase/test5.vcf
echo "PASSED"

echo "Test 6"
./akt pedphase test/pedphase/test6.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test6.out
diff test/pedphase/test6.out test/pedphase/test6.vcf
./akt pedphase test/pedphase/test6.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test6.out
diff test/pedphase/test6.out test/pedphase/test6.vcf
echo "PASSED"

echo "Test 7"
./akt pedphase test/pedphase/test7.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test7.out
diff test/pedphase/test7.out test/pedphase/test7.vcf
./akt pedphase test/pedphase/test7.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test7.out
diff test/pedphase/test7.out test/pedphase/test7.vcf
echo "PASSED"

echo "Test 8"
./akt pedphase test/pedphase/test8.vcf.gz -p test/pedphase/pedigree.fam -o - 2> /dev/null | bcftools view -H > test/pedphase/test8.out
diff test/pedphase/test8.out test/pedphase/test8.vcf
./akt pedphase test/pedphase/test8.vcf.gz -o - 2> /dev/null | bcftools view -H > test/pedphase/test8.out
diff test/pedphase/test8.out test/pedphase/test8.vcf
echo "PASSED"

echo "ALL TESTS PASSED!"
