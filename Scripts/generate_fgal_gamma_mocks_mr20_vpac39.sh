#! /bin/bash

output=stomp_sdssmock_main20_esmeralda
maxtheta=0.1
binning=10
halo=/ssd1/Research/halo_files/Esmeralda
fof=/home/piscioja/halobias/halobias_fof_nfw
so=/home/piscioja/halobias/halobias_so_nfw

map_dir=/home/piscioja/SDSSPix/Maps
mapfile0=$map_dir/window.dr72brvoid0.stripe_trim.pix
mapfile1=$map_dir/lss_geometry.dr72.stripe_trim.pix
mapfile3=$map_dir/lss_geometry.fullsphere.pix
map_overlap2=$map_dir/window.dr72brvoid1.overlap.noweight.pix 


stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=11.8013669053163
siglogM=0.0973817154668221
logM0=12.687852882335
logM1=13.0400932413139
alpha=1.05768884410071
delta_vir=200.0 #linking length of 0.2
Mstar=2.29E12
redshift=0.082

padtowidth=3



echo "$1"

i=$1

mkdir /ssd1/Research/GalaxyCatalogues/Esmeralda/$i
cd /ssd1/Research/GalaxyCatalogues/Esmeralda/$i

for gamma1 in $(seq 1.9 0.008 2.9 )
do



for fgal1 in $(seq 0.21 0.005 0.41)
do



gamma="$(printf "%5.4f" "$gamma1")"
fgal="$(printf "%5.4f" "$fgal1")"

echo "$gamma $fgal"


time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < $halo/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > $halo/fff/halobias_fof_nfw_${i}_fff 

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.106 $mapfile1 0.6 0 < $halo/fff/halobias_fof_nfw_${i}_fff > halobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.f$fgal.g$gamma.fof.lss_geometry  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 


rm $halo/fff/halobias_fof_nfw_${i}_fff
rm halobias_fof_nfw_${i}_galaxies




done
done




