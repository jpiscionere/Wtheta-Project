#! /bin/bash

output=stomp_sdssmock_main20_esmeralda
maxtheta=0.1
binning=10
fof=/home/piscioja/halobias/halobias_fof_nfw
so=/home/piscioja/halobias/halobias_so_nfw

map_dir=/data2/jap/SDSSPix/Maps
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


for i in $(seq 3024 1 3030)
do
cd /data2/jap/GalaxyCatalogues/FgalGamma/Esmeralda/$i

for gamma1 in $(seq 1.5 0.02 2.5 )
do



for fgal1 in $(seq 0.01 0.01 0.51)
do



gamma="$(printf "%4.3f" "$gamma1")"
fgal="$(printf "%4.3f" "$fgal1")"

echo "$gamma $fgal"


time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < /data0/LasDamas/Esmeralda/$i/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > ehalobias_fof_nfw_${i}_fff 
#time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < /data2/jap/GalaxyCatalogues/test> ehalobias_fof_nfw${i}_fff 

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.106 /data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix 0.6 0 < ehalobias_fof_nfw_${i}_fff > ehalobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=ehalobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.f$fgal.g$gamma.fof.lss_geometry  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 

rm ehalobias*




done
done
done



