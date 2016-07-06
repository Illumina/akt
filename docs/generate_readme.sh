#!/bin/bash


echo "#Using akt" >> usage.md
cat Using.md >> usage.md; 		../akt 2> tmp; 				sed -e '/replace_akt_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat kin.md >> usage.md; 		../akt kin 2> tmp; 			sed -e '/replace_kin_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat relatives.md >> usage.md; 	../akt relatives 2> tmp; 	sed -e '/replace_relatives_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat mendel.md >> usage.md; 	../akt mendel 2> tmp; 	sed -e '/replace_mendel_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat pca.md >> usage.md; 		../akt pca 2> tmp; 			sed -e '/replace_pca_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat cluster.md >> usage.md; 	../akt cluster 2> tmp; 		sed -e '/replace_cluster_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat stats.md >> usage.md; 		../akt stats 2> tmp; 		sed -e '/replace_stats_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat ibd.md >> usage.md; 		../akt ibd 2> tmp; 			sed -e '/replace_ibd_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat LDplot.md >> usage.md; 		../akt LDplot 2> tmp; 		sed -e '/replace_LDplot_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat admix.md >> usage.md; 		../akt admix 2> tmp; 		sed -e '/replace_admix_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 
cat metafreq.md >> usage.md; 		../akt metafreq 2> tmp; 		sed -e '/replace_metafreq_run/{' -e 'r tmp' -e 'd' -e '}' -i usage.md; 

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


declare -A options;
options=(["-o"]="output_option"); 
options+=(["-O"]="outputfmt_option");
options+=(["-r"]="regions_option");
options+=(["-R"]="regionsfile_option");
options+=(["-s"]="samples_option");
options+=(["-S"]="samplesfile_option");
options+=(["-T"]="targetfile_option");
options+=(["-t"]="target_option");
options+=(["-n"]="nthread_option");
options+=(["-a"]="aftag_option");
options+=(["-f"]="pairfile_option");
options+=(["-h"]="thin_option");
options+=(["-m"]="maf_option");

for i in "${!options[@]}"; do
   line=`grep "^$i" common_options.md`
   sed -i -e 's/'"${options["$i"]}"'/'"$line"'/g' examples.md
   sed -i -e 's/'"${options["$i"]}"'/'"$line"'/g' usage.md
done

rm tmp;
