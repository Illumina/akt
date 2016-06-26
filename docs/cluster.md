##cluster
```
replace_cluster_run
```
* -k : Number of clusters. Examine the data to guess this or analyse silhouette scores.
* -i : Random seed.
* -a : 0 to use k++ means clustering, 1 to use EM with Gaussians, 2 to use density method.
* -C : file with initial guess for cluster centres. If analysing 2d data this should contain 2 columns and K rows.
* -c : Which columns in input file to use.
* -I : Maximum number of iterations to use for alg 0 or 1.
* -d : radius around each point for counting density.
* -p : min density for cluster centre.
* -D : min radius for cluster centre.
* --density-plot : plot delta-density graph and finish.
* -e : Calculate silhouette score (goodness of cluster assignment).

e.g. cluster the first three principle components from the output of vcfpca using k++ means. 
Different clusters should correspond to different ancestries
```
./akt cluster projections -k 4 -c 2-4 -e > clustered
```
Each block in `clustered` contains all the data for points in that cluster with the row layout
```
P2 P3 P4 (SILHOUETTE) CLUSTERID P1 P5 ...
```
To visually inspect the clustered data try the following gnuplot command
```
gnuplot> splot for [i=0:3:1] 'clustered' every:::i::i
```
The silhouette column can be used to filter on badly clustered samples. 
~1 means a sample is closer to the centroid of its assigned cluster than the next nearest cluster.
A sample can be badly clustered for many reasons: bad choice of K, convergence of the clustering algorithm 
to bad local minimum, individual is of mixed ancestry. You could filter out badly clustered samples using
```
awk '{ if($4 > 0.5) print $0}' clustered > well_clustered
```
###density clustering
This is a better clustering algorithm but requires some tuning by hand. First run
```
./akt cluster projections -c 2-4 -a 2 -d 1 --density-plot > dplot
```
the file dplot is a two column data file of density versus delta 
(the distance to the closest data point of higher density). By examining this plot outliers 
are obvious and a cutoff in density and delta can be determined that isolates these points. 
Once this cutoff is determined run
```
./akt cluster projections -c 2-4 -a 2 -d 1 -p 0 -D 75 > clustered
```
Where 0 and 75 are the cutoffs determined from `dplot`. This will output the number of clusters that were found 
(including one for unassigned points) and `clustered` will contain the clustered data in the usual format.
```
gnuplot> splot for [i=0:4:1] 'clustered' every:::i::i
```
Unassigned data points are always put in cluster 0. The silhouette score for unassigned data is undefined.
