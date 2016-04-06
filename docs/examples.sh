#!/bin/bash
if [ ! -f ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf ]
then
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf
wget https://s3-eu-west-1.amazonaws.com/akt-examples/1000G/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi
fi

##Discovering Cryptic Relations
../akt kin ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R ../data/1000G.snps.nochr.vcf.gz -n 4 > test_kin
../akt relatives test_kin -p test -g > test_relatives
neato test.allgraph -Tpng -O
dot test.Fam0.graph -Tpng -O
dot test.Fam134.graph -Tpng -O
grep '^Unrel' test_relatives | awk '{print $2}' > test.unrel

##Discovering Ancestry
../akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R ../data/1000G.snps.nochr.vcf.gz > test_pca
../akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -R ../data/1000G.snps.nochr.vcf.gz -S phase1.unrel > test_pca.unrel
../akt pca ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -w ../data/1000G.snps.nochr.vcf.gz > test_pcaproj
../akt cluster test_pca -c 2-3 -a 2 -d 1 --density-plot > test_dplot
../akt cluster test_pca -c 2-3 -a 2 -d 1 -p -1 -D 700 > test_allclusters
../akt cluster test_pca -c 2-3 -a 2 -d 1 -p 0 -D 700 > test_allclusters1

##Allele Frequency and Correlation
grep Cluster1 test_allclusters | awk '{print $4}' > test.Cluster1
../akt stats ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -a "C1" -O b -o ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.sitesx.bcf -S test.Cluster1 -r 20 -x
../akt stats ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -a "EAS" -O b -o sites.EAS.bcf -S test.Cluster1 -r 20 -F 10000
bcftools index sites.EAS.bcf


##Segment Sharing
echo HG00464 HG00537 > pairs.txt
echo HG00477 HG00544 >> pairs.txt
echo HG00477 HG00542 >> pairs.txt
echo HG00475 HG00544 >> pairs.txt
echo HG00475 HG00542 >> pairs.txt
../akt ibd ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -T sites.EAS.bcf -a "EAS" -f pairs.txt > test_ibd

#BEAGLE, really slow!
#bcftools view ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -S test.Cluster1 -O z -o test.20.vcf.gz --force-samples;
#java -jar ibdseq.r1206.jar gt=test.20.vcf.gz out=Seq.EAS.test.20 nthreads=4 

##Correlation Matrix
../akt LDplot ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20:0-200000 > test_sigma

##Calculating Admixture Fractions
echo 40 -4 AFR 1 0 > centre.txt
echo -12 21 EUR 0 1 >> centre.txt
echo -19 -27 EAS 0 0 >> centre.txt
../akt admix test_pcaproj -c 2-3 -C centre.txt > test_admix
../akt admix test_pcaproj -c 2-6 -C ../data/1000G.pca_to_admix > test_alladmix

##Make plots
##NB: docs/AMR.samples does not equal data/AMR.samples!
for pop in AMR AFR EAS EUR SAS; 
do 
	while read line; do grep "$line" test_pca; done < "$pop".samples > test_pca_"$pop"; 

done
for pop in AMR AFR EAS EUR; 
do 
	while read line; do grep "$line" test_pcaproj; done < "$pop".samples > test_pcaproj_"$pop"; 
	while read line; do grep "$line" test_admix; done < "$pop".samples > test_admix_"$pop"; 
	while read line; do grep "$line" test_alladmix; done < "$pop".samples > test_alladmix_"$pop"; 

done
gnuplot -e "load 'plot.gnu'"

##Makes the complicated ibd segment plot
cp ibd.gnu tmp.gnu
i=1;
l=0;
while read line; 
do 
	p1=`echo $line | awk '{print $1}'`; 
	p2=`echo $line | awk '{print $2}'`;
printf "#%s %s\n" $p1, $p2 >> tmp.gnu
printf "set label '%s, %s' at 95000000,%s.1 right\n" $p1 $p2 $l >> tmp.gnu

grep $p1 test_ibd | grep $p2 | awk '{ printf "set arrow %s from %s,'"$l"' to %s,'"$l"' nohead\n", '"$i"' + a, $3, $4; a++ }' >> tmp.gnu
grep $p1 Seq.EAS.test.20.ibd | grep $p2 | awk '{ if($7-$6>1000000){ printf "set arrow %s from %s,'"$l"'.2 to %s,'"$l"'.2 nohead ls 1\n", '"$i"'+50+a, $6, $7; a++} }' >> tmp.gnu

l=$(($l+1))
i=$(($i+100))
done < pairs.txt

printf "set term pngcairo\nset out \'segment.png\'\n" >> tmp.gnu

echo plot -10 notitle >> tmp.gnu
printf "#EOF\n" >> tmp.gnu
gnuplot -e "load 'tmp.gnu'"
rm tmp.gnu
