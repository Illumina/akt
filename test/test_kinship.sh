reg=../data/wgs.grch37.vcf.gz
data=ALL.cgi_multi_sample.20130725.pruned.snps.bcf

##calculate kinship coefficients
time ../akt kin -M 0 -@ 4 -R $reg $data > kinship0.txt
time ../akt kin -M 1 -@ 4 -R $reg $data > kinship1.txt
time ../akt kin -@ 4 -F $reg $data > kinship.txt

# check for unrelated
../akt unrelated kinship1.txt | sort > unrelated.out
#diff unrelated.ids unrelated.out #fixme: make this deterministic

##find relatives
time ../akt relatives -p n433 kinship.txt > relatives.out
python ped_compare.py  n433.fam  20130606_g1k.fam
