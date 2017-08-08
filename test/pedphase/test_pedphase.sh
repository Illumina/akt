#!/usr/bin/env bash

echo "Test 1"
../akt pedphase pedphase/test1.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test1.out
diff pedphase/test1.out pedphase/test1.vcf
../akt pedphase pedphase/test1.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test1.out
diff pedphase/test1.out pedphase/test1.vcf
echo "PASSED"

echo "Test 2"
../akt pedphase pedphase/test2.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test2.out
diff pedphase/test2.out pedphase/test2.vcf
../akt pedphase pedphase/test2.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test2.out
diff pedphase/test2.out pedphase/test2.vcf
echo "PASSED"

echo "Test 3"
../akt pedphase pedphase/test3.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test3.out
diff pedphase/test3.out pedphase/test3.vcf
../akt pedphase pedphase/test3.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test3.out
diff pedphase/test3.out pedphase/test3.vcf
echo "PASSED"

echo "Test 5"
../akt pedphase pedphase/test5.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test5.out
diff pedphase/test5.out pedphase/test5.vcf
../akt pedphase pedphase/test5.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test5.out
diff pedphase/test5.out pedphase/test5.vcf
echo "PASSED"

echo "Test 6"
../akt pedphase pedphase/test6.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test6.out
diff pedphase/test6.out pedphase/test6.vcf
../akt pedphase pedphase/test6.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test6.out
diff pedphase/test6.out pedphase/test6.vcf
echo "PASSED"

echo "Test 7"
../akt pedphase pedphase/test7.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test7.out
diff pedphase/test7.out pedphase/test7.vcf
../akt pedphase pedphase/test7.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test7.out
diff pedphase/test7.out pedphase/test7.vcf
echo "PASSED"

echo "Test 8"
../akt pedphase pedphase/test8.vcf.gz -p pedphase/pedigree.fam -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test8.out
diff pedphase/test8.out pedphase/test8.vcf
../akt pedphase pedphase/test8.vcf.gz -o - 2> /dev/null | grep -A1000 CHROM > pedphase/test8.out
diff pedphase/test8.out pedphase/test8.vcf
echo "PASSED"

echo "ALL TESTS PASSED!"
