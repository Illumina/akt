#!/bin/bash

echo -e "# akt - ancestry and kinship toolkit" > ../README.md
cat License.md >> ../README.md
cat Install.md >> ../README.md
echo >> ../README.md
echo "##Using akt" >> ../README.md
cat Using.md >> ../README.md; 		../akt 2> tmp; 				sed -e '/replace_akt_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat kin.md >> ../README.md; 		../akt kin 2> tmp; 			sed -e '/replace_kin_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat relatives.md >> ../README.md; 	../akt relatives 2> tmp; 	sed -e '/replace_relatives_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat mendel.md >> ../README.md; 	../akt mendel 2> tmp; 	sed -e '/replace_mendel_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat pca.md >> ../README.md; 		../akt pca 2> tmp; 			sed -e '/replace_pca_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat cluster.md >> ../README.md; 	../akt cluster 2> tmp; 		sed -e '/replace_cluster_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat stats.md >> ../README.md; 		../akt stats 2> tmp; 		sed -e '/replace_stats_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat ibd.md >> ../README.md; 		../akt ibd 2> tmp; 			sed -e '/replace_ibd_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat LDplot.md >> ../README.md; 		../akt LDplot 2> tmp; 		sed -e '/replace_LDplot_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 
cat admix.md >> ../README.md; 		../akt admix 2> tmp; 		sed -e '/replace_admix_run/{' -e 'r tmp' -e 'd' -e '}' -i ../README.md; 

echo -e "##Example Workflow" >> ../README.md
cat Data.md >> ../README.md
cat Cryptic.md >> ../README.md
cat Mendel_example.md >> ../README.md
cat Ancestry.md >> ../README.md
cat AF.md >> ../README.md
cat SegIBD.md >> ../README.md
cat Corplot.md >> ../README.md
cat AdFrac.md >> ../README.md

declare -A options;
options=(["-o"]="output_option"); 
options+=(["-O"]="outputfmt_option");
options+=(["-r"]="regions_option");
options+=(["-R"]="regionsfile_option");
options+=(["-s"]="samples_option");
options+=(["-S"]="samplesfile_option");
options+=(["-T"]="targetfile_option");
options+=(["-n"]="nthread_option");
options+=(["-a"]="aftag_option");
options+=(["-f"]="pairfile_option");
options+=(["-h"]="thin_option");
options+=(["-m"]="maf_option");

for i in "${!options[@]}"; do
   line=`grep "^$i" common_options.md`
   sed -i -e 's/'"${options["$i"]}"'/'"$line"'/g' ../README.md
done

rm tmp;
