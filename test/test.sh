#this is a basic set of tests for akt.

##build akt
cd ../
make clean
make -j 4
cd test/

##get data
data=ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
if [ ! -f $data ]
then
    wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/${data}
fi

if [ ! -f ${data}.csi ]
then
    wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/${data}.csi
fi

##pca of data
time ../akt pca -R ../data/1000G.snps.nochr.vcf.gz $data  > pca1.txt

##project data onto 1000G PCs
time ../akt pca -w ../data/1000G.snps.nochr.vcf.gz $data  > pca2.txt
Rscript ../scripts/1000G_pca.R pca2.txt 

##calculate kinship coefficients
time ../akt kin -n 4 -R ../data/1000G.snps.nochr.vcf.gz $data > kinship.txt

##find relatives
time ../akt relatives -p n433 kinship.txt > relatives.out

##check mendel error rates of relatives
time ../akt mendel -p n433.fam $data > mendel.out
Rscript  ../scripts/mendel.R mendel.out

##clustering on pca plot
time ../akt cluster pca1.txt -c 2-3 -a 2 -d 1 --density-plot > dplot.out
time ../akt cluster pca1.txt -c 2-3 -a 2 -d 1 -p -1 -D 700 > cluster.out

##calculate stats for chrom 20
time ../akt stats $data -O b -o sites.bcf -r 20
bcftools index sites.bcf

##calculate ibd sharing
time ../akt ibd $data -r 20 -T sites.bcf -s "HG00475,HG00542" > ibd.out

##create LD plot
time ../akt LDplot $data -r 20:0-200000 > LDplot.out

##calculate admixture fractions
time ../akt admix pca2.txt -c 2-6 -C ../data/1000G.pca_to_admix > admix.out
