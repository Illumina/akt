The site-only vcfs in this directory contain a sparse set (one per ~200kb) of 17535 "reliable" SNPs.  These are reasonably common and consistently called by different bioinformatics pipelines. They should be useful for PCA and relatedness analyses which require LD-pruned variants. We have one site list appropriate for whole-genome sequence (WGS) data and one approprate for whole-exome sequence data (WES).

There is one vcf where chromosomes have the chr prefix (hg19/hg38) and one vcf where the prefix is removed (grch37/grch38). 

**WGS data:**
* bi-allelic SNPs in 1000G Phase 3 with global MAF>=5%
* present on the Human Core Exome Microarray
* present on the Infinium Core Exome Microarray
* present on the omni2.5 microarray
* present on the affy6 microarray
* present on v4 of 23andMe microarray
* present in high-coverage samples sequenced at Illumina with similar allele frequency to 1000G

**WES data:**
* bi-allelic SNPs in 1000G Phase 3 with global MAF>=5%
* present in ExAC
* present in high-coverage samples sequenced at Illumina with similar allele frequency to 1000G
* LD pruned 



The INFO also contains the PCA coefficients for 1000G Phase3 which allows new samples to be easily projected onto the 1000G principal components. This is only applicable for the grch37/hg19 site lists, since 1000G is not available in b38.

The files `wgs.1000G.phase3.pca` and `wes.1000G.phase3.pca` contain the first twenty principal components calculated on the 2504 samples in the 1000G Phase 3 release using these SNPs.

You can generate a simple plot in R of the 1000G WGS PCs via:

```
mycol <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
ogp.sample <- read.table("integrated_call_samples_v3.20130502.ALL.panel",header=T,as.is=T)
ogp.pc <- read.table("wgs.1000G.phase3.pca",as.is=TRUE)
ogp.pc$pop <- ogp.sample$super_pop[match(ogp.sample[,1],ogp.pc[,1])]
names(mycol) <- unique(ogp.sample$super_pop)

par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(ogp.pc[,2]~ogp.pc[,3],col=mycol[ogp.sample$super_pop],xlab='PC2',ylab='PC1')
plot(ogp.pc[,4]~ogp.pc[,3],col=mycol[ogp.sample$super_pop],xlab='PC2',ylab='PC3')
```

