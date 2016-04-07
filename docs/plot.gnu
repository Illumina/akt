#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 4    last modified 2013-10-02 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal wxt 0
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 25, 64, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ * : * ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_GB.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
GNUTERM = "wxt"
set term pngcairo; 


set out 'test_kin.png'
set xlabel 'IBD0'
set ylabel 'IBD1'
plot 'test_kin' using 3:4 ti 'pair of samples ='


set out 'test_pca12.png'; 
set xlabel 'PCA1'
set ylabel 'PCA2'
plot 'test_pca_SAS' using 2:3 ti 'SAS',\
'test_pca_EUR' using 2:3 ti 'EUR',\
'test_pca_EAS' using 2:3 ti 'EAS',\
'test_pca_AMR' using 2:3 ti 'AMR',\
'test_pca_AFR' using 2:3 ti 'AFR'

set out 'test_pca32.png'; 
set xlabel 'PCA3'
set ylabel 'PCA2'
plot 'test_pca_SAS' using 4:3 ti 'SAS',\
'test_pca_EUR' using 4:3 ti 'EUR',\
'test_pca_EAS' using 4:3 ti 'EAS',\
'test_pca_AMR' using 4:3 ti 'AMR',\
'test_pca_AFR' using 4:3 ti 'AFR'

set out 'test_cluster.png'; 
set xlabel 'PCA1'
set ylabel 'PCA2'
plot 'test_allclusters' using 1:2 every :::0::0 ti 'Cluster1' lt 3,\
'test_allclusters' using 1:2 every :::1::1 ti 'Cluster2' lt 4,\
'test_allclusters' using 1:2 every :::2::2 ti 'Cluster3' lt 1,\
'test_allclusters' using 1:2 every :::3::3 ti 'Cluster4' lt 2,\
'test_allclusters' using 1:2 every :::4::4 ti 'Cluster5' lt 5

set out 'test_cluster1.png'; 
set xlabel 'PCA1'
set ylabel 'PCA2'
plot 'test_allclusters1' using 1:2 every :::0::0 ti 'Cluster0' lt 6,\
'test_allclusters1' using 1:2 every :::1::1 ti 'Cluster1' lt 3,\
'test_allclusters1' using 1:2 every :::2::2 ti 'Cluster2' lt 4,\
'test_allclusters1' using 1:2 every :::3::3 ti 'Cluster3' lt 1,\
'test_allclusters1' using 1:2 every :::4::4 ti 'Cluster4' lt 2,\
'test_allclusters1' using 1:2 every :::5::5 ti 'Cluster5' lt 5

set out 'test_pcaproj12.png'; 
set xlabel 'PCA1'
set ylabel 'PCA2'
plot 'test_pcaproj_EUR' using 2:3 ti 'EUR',\
'test_pcaproj_EAS' using 2:3 ti 'EAS',\
'test_pcaproj_AMR' using 2:3 ti 'AMR',\
'test_pcaproj_AFR' using 2:3 ti 'AFR'

set out 'test_admix.png'; 
set xlabel '%AFR'
set ylabel '%EUR'
plot 'test_admix_EUR' using 2:3 ti 'EUR',\
'test_admix_EAS' using 2:3 ti 'EAS',\
'test_admix_AMR' using 2:3 ti 'AMR',\
'test_admix_AFR' using 2:3 ti 'AFR'

set out 'test_sigma.png'
set xrange [800:1200]
set yrange [800:1200] reverse
plot 'test_sigma' matrix with image
unset yrange
set auto

set out 'test_alladmix.png'; 
set multiplot layout 2,2 rowsfirst; 
unset key; set xlabel '%EUR'; set ylabel '%AMR'; plot 'test_alladmix_AFR' using 2:3 ti 'AFR', 'test_alladmix_AMR' using 2:3 ti 'AMR', 'test_alladmix_EUR' using 2:3 ti 'EUR', 'test_alladmix_EAS' using 2:3 ti 'EAS';
unset key; set xlabel '%EUR'; set ylabel '%AFR'; plot 'test_alladmix_AFR' using 2:4 ti 'AFR', 'test_alladmix_AMR' using 2:4 ti 'AMR', 'test_alladmix_EUR' using 2:4 ti 'EUR', 'test_alladmix_EAS' using 2:4 ti 'EAS';
unset key; set xlabel '%EUR'; set ylabel '%EAS'; plot 'test_alladmix_AFR' using 2:5 ti 'AFR', 'test_alladmix_AMR' using 2:5 ti 'AMR', 'test_alladmix_EUR' using 2:5 ti 'EUR', 'test_alladmix_EAS' using 2:5 ti 'EAS';
unset key; set xlabel '%EUR'; set ylabel '%SAS'; plot 'test_alladmix_AFR' using 2:6 ti 'AFR', 'test_alladmix_AMR' using 2:6 ti 'AMR', 'test_alladmix_EUR' using 2:6 ti 'EUR', 'test_alladmix_EAS' using 2:6 ti 'EAS'


#    EOF
