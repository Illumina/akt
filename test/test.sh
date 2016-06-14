#this is a basic set of tests for akt.


reg=../data/wgs.grch37.vcf.gz
##get data
data=ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/${data}
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/${data}.csi

##pca of data
time ../akt pca -R $reg $data  > pca1.txt

##project data onto 1000G PCs
time ../akt pca -w $reg $data  > pca2.txt
Rscript ../scripts/1000G_pca.R pca2.txt 

##calculate kinship coefficients
time ../akt kin -M 0 -n 4 -R $reg $data > kinship0.txt
time ../akt kin -M 1 -n 4 -R $reg $data > kinship1.txt
time ../akt kin -M 2 -n 4 -R $reg $data > kinship2.txt
time ../akt kin -n 4 -F $reg $data > kinship.txt


##find relatives
time ../akt relatives -p n433 kinship.txt > relatives.out
python ped_compare.py  n433.fam  20130606_g1k.fam

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

##simple metafreq test
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/UK10K_COHORT.chr20.20140722.sites.vcf.gz
wget --continue https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/UK10K_COHORT.chr20.20140722.sites.vcf.gz.tbi

time ../akt metafreq UK10K_COHORT.chr20.20140722.sites.vcf.gz ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -Oz -o uk10_1000g.frq.chr20.vcf.gz
