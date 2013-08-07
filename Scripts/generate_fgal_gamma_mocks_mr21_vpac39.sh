#! /bin/bash

output=stomp_sdssmock_main21_carmen
maxtheta=.10
binning=10
halo=/ssd1/Research/halo_files/Carmen
fof=/home/piscioja/halobias/halobias_fof_nfw
so=/home/piscioja/halobias/halobias_so_nfw


map_dir=/home/piscioja/SDSSPix/Maps
mapfile0=$map_dir/window.dr72brvoid0.stripe_trim.pix
mapfile1=$map_dir/lss_geometry.dr72.stripe_trim.pix
mapfile3=$map_dir/lss_geometry.fullsphere.pix
map_overlap2=$map_dir/window.dr72brvoid1.overlap.noweight.pix

stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=12.76422334
siglogM=0.71194502
logM0=13.07148775
logM1=13.70323156
alpha=1.07455104
delta_vir=200.0 #linking length of 0.2
Mstar=1.97E12
redshift=0.132


i=$1
echo "$1"




mkdir /ssd1/Research/GalaxyCatalogues/Carmen/$i
cd /ssd1/Research/GalaxyCatalogues/Carmen/$i

for gamma1 in $(seq 1.8 0.006 2.1 )
do


for fgal1 in $(seq 0.01 0.0025 0.11)
do


gamma="$(printf "%5.4f" "$gamma1")"
fgal="$(printf "%5.4f" "$fgal1")"

echo "$gamma $fgal"


time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -21 $delta_vir $Mstar $redshift  trashfile.out -$i < $halo/Carmen_${i}_z0p132_fof_b0p2.dpp.halos > $halo/fff/halobias_fof_nfw_${i}_fff 

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.165 $mapfile1 0.6 0 < $halo/fff/halobias_fof_nfw_${i}_fff > halobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.f$fgal.g$gamma.fof.lss_geometry  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 

rm $halo/fff/halobias_fof_nfw_${i}_fff
rm halobias*

done
done






