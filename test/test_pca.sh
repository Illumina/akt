reg=../data/wgs.grch37.vcf.gz
data=ALL.cgi_multi_sample.20130725.pruned.snps.bcf
##pca of data
time ../akt pca -R $reg $data  > pca1.txt

##project data onto 1000G PCs
time ../akt pca -W $reg $data  > pca2.txt
Rscript ../scripts/1000G_pca.R pca2.txt 
