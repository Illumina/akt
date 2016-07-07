#!/bin/bash
rm usage.md
cat Using.md >> usage.md; 		../akt 2> tmp; 				sed -e '/replace_akt_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat pca.md >> usage.md; 		../akt pca 2> tmp; 			sed -e '/replace_pca_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat kin.md >> usage.md; 		../akt kin 2> tmp; 			sed -e '/replace_kin_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat relatives.md >> usage.md; 	../akt relatives 2> tmp; 	sed -e '/replace_relatives_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat mendel.md >> usage.md; 	../akt mendel 2> tmp; 	sed -e '/replace_mendel_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat cluster.md >> usage.md; 	../akt cluster 2> tmp; 		sed -e '/replace_cluster_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat stats.md >> usage.md; 		../akt stats 2> tmp; 		sed -e '/replace_stats_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat ibd.md >> usage.md; 		../akt ibd 2> tmp; 			sed -e '/replace_ibd_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat LDplot.md >> usage.md; 		../akt LDplot 2> tmp; 		sed -e '/replace_LDplot_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat admix.md >> usage.md; 		../akt admix 2> tmp; 		sed -e '/replace_admix_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat metafreq.md >> usage.md; 		../akt metafreq 2> tmp; 		sed -e '/replace_metafreq_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 


sed -i s,"output_option","**-o** *FILE* see common options  ",g usage.md
sed -i s,"outputfmt_option","**-O** *z|b|v|u* see common options  ",g usage.md  
sed -i s,"regionsfile_option","**-R** *FILE* see common options  ",g usage.md  
sed -i s,"regions_option","**-r** *REGION* see common options  ",g usage.md  
sed -i s,"targetfile_option","**-T** *FILE* see common options  ",g usage.md  
sed -i s,"target_option","**-t** *TARGET* see common options  ",g usage.md  
sed -i s,"samplesfile_option","**-S** *SAMPLES* see common options  ",g usage.md 
sed -i s,"samples_option","**-s** *SAMPLES* see common options  ",g usage.md   
sed -i s,"nthread_option","**-n** *VALUE* see common options  ",g usage.md  
sed -i s,"maf_option","**-m** *VALUE* see common options  ",g usage.md

rm tmp;


echo -e "#Example Workflow" >> examples.md
cat Data.md >> examples.md
cat Cryptic.md >> examples.md
cat Mendel_example.md >> examples.md
cat Ancestry.md >> examples.md
cat AF.md >> examples.md
cat SegIBD.md >> examples.md
cat Corplot.md >> examples.md
cat AdFrac.md >> examples.md
cat MF.md >> examples.md

