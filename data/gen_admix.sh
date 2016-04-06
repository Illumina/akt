#!/bin/bash

##All of phase3 1000Genomes
data=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf
if [ ! -f $data ]
then
	echo "DOWNLOADING ALL 1000G phase3 vcfs - this will take a while!"
	files=""
	for((i=1;i<=22;i++));
	do
		files="$files ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
		wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
	done
	#bcftools concat $files -O b -o ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf
fi
if [ ! -f ${data}.csi ]
then
    bcftools index $data
fi

##Project phase3 onto principle components in data
../akt pca $data -w 1000G.snps.nochr.vcf.gz > phase3_proj_pca

##use known ancestry to fill in admixture transformation
#Note: AMR.samples is a subset of 1000G AMR individuals
i=0
for pop in EUR AMR AFR EAS SAS; 
do 
	while read line; 
	do 
		grep "$line" phase3_proj_pca
	done < "$pop".samples > phase3_proj_pca_"$pop";

	../akt cluster -k 1 -o centre_"$pop" -c 2-6 phase3_proj_pca_"$pop" &> junk 
	
	rm junk
	rm centre_"$pop"
	rm phase3_proj_pca_"$pop";
	
	awk '{printf "%s\t%s\t", $0, "'"$pop"'"; for(j=0; j<NF; j++){ if(j=="'"$i"'"){printf "1\t";}else{printf "0\t";} } printf "\n" }' centre_"$pop"
	i=$(($i+1))
	
done > 1000G.pca_to_admix

rm phase3_proj_pca
