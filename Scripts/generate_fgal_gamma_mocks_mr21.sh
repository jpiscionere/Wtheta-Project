#! /bin/bash

output=stomp_sdssmock_main21_carmen
maxtheta=.10
binning=10
fof=/home/piscioja/halobias/halobias_fof_nfw
so=/home/piscioja/halobias/halobias_so_nfw
mapfile0=/home/piscioja/GalaxyCatalogues/window.dr72brvoid0.stripe_trim.pix
mapfile1=/data2/cameron/for_jen/Archive/lss_geometry.dr72.stripe_trim.pix
mapfile3=/data2/jap/GalaxyCatalogues/lss_geometry.fullsphere.pix
mapfile4=/data2/jap/GalaxyCatalogues/lss_geometry.fullsphere_single_index.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=12.76422334
siglogM=0.71194502
logM0=13.07148775
logM1=13.70323156
alpha=1.07455104
delta_vir=200.0 #linking length of 0.2
Mstar=1.97E12
redshift=0.132





for i in $(seq 2024 1 2039)
do
cd /data2/jap/GalaxyCatalogues/FgalGamma/Carmen/$i

for gamma1 in $(seq 1.5 0.02 2.5 )
do


for fgal1 in $(seq 0.01 0.01 0.51)
do


gamma="$(printf "%4.3f" "$gamma1")"
fgal="$(printf "%4.3f" "$fgal1")"

echo "$gamma $fgal"


time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -21 $delta_vir $Mstar $redshift  trashfile.out -$i < /data0/LasDamas/Carmen/$i/Carmen_${i}_z0p132_fof_b0p2.dpp.halos > hhalobias_fof_nfw_${i}_fff 

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.165 /data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix 0.6 0 < hhalobias_fof_nfw_${i}_fff > hhalobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=hhalobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.f$fgal.g$gamma.fof.lss_geometry  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 
#time $stompdirectory/stomp_galaxy_autocorrelation_jack --map_file=$mapfile0 --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.f$fgal.g$gamma.fof.lss_geometry.jack  --theta_max=$maxtheta --n_bins_per_decade=$binning n_random=1 --single_index

rm hhalobias*

done
done
done





