##Segment sharing
We can try to find long regions shared IBD between two samples. 
We choose to demonstrate this for samples which are closely related, though not flagged by the kin tool or by 1000 Genomes. 
These are all East Asian individuals. 
```
cat pairs.txt
HG00464 HG00537
HG00477 HG00544
HG00477 HG00542
HG00475 HG00544
HG00475 HG00542
```
Then run
```
./akt ibd ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -T sites.EAS.bcf -a "EAS" -f pairs.txt > test_ibd
```
We only look at chromosome 20 and a small number of pairs. `akt ibd` uses the LD score metric computed by `akt stats`. If
this is not present the markers are assumed to be in linkage equilibrium. Unfortunately there isn't good truth data for IBD. 
We run  [IBDseq](http://faculty.washington.edu/browning/ibdseq.html) on the same data
```
bcftools view ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -r 20 -S test.Cluster1 -O z -o test.20.vcf.gz;
java -jar ibdseq.r1206.jar gt=EAS.test.20.vcf.gz out=Seq.EAS.test.20
```
Comparing the IBD segments gives

![alt text](https://git.illumina.com/Bioinformatics/akt/blob/master/docs/segment.png)

Where the yscale is arbitrary and the xscale is position along chromosome 20. We show the segments 
as reported by IBDseq in red and by akt in black. For each of the 5 pairs examined 
(chosen because of the presence of significant IBD sharing) the two methods give similar answers. 
We expect IBDseq to be more accurate since it uses phased haplotypes where akt uses only genotypes.
Compared to IBDseq, akt produces some false positives and some broken segments reported as continuous 
by IBDSeq. These can be controlled by tuning the score and minimum segment length thresholds.
