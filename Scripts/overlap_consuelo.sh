#! /bin/bash
fof=/home/piscioja/halobias/halobias_fof_nfw
output=stomp_sdssmock_main18_consuelo
maxtheta=10
binning=10
mapfile1=/data2/jap/GalaxyCatalogues/lss_geometry.dr72.stripe_trim.pix
mapfile3=/data2/jap/GalaxyCatalogues/lss_geometry.fullsphere.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=11.18
siglogM=0.19
logM0=9.81
logM1=12.42
alpha=1.04
delta_vir=377.0 #linking length of 0.2
Mstar=2.49E12
redshift=0.054
gamma=1
fgal=1


#$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=/data2/jap/GalaxyCatalogues/Randoms/sdssmock_gamma_main21.rand_10x_ns.rdcz.dat  -output_tag=sdssmock_gamma_main21.RR  --theta_max=$maxtheta --n_bins_per_decade=$binning


for i in $(seq 4005 1 4015)
do


time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -19.0 $delta_vir $Mstar $redshift trashfile.out -$i < /data0/LasDamas/Consuelo/$i/Consuelo_${i}_z0p054_fof_b0p156.fdpp.halos > halobias_fof_nfw_${i}_fff_18


time /home/piscioja/halobias/makemock 0 1 0 0 0 0 0.02 0.042 /data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix 0.6 0 < halobias_fof_nfw_${i}_fff_18 > halobias_nfw_nopart_${i}_galaxies_18





for overlap in n_o s_o d_o n_c


do

if	[[ $overlap == "n_o" ]]
then
	./IDfib2 halobias_nfw_nopart_${i}_galaxies_18 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_${overlap}_18

elif	[[ $overlap == "s_o" ]]
then
	./IDfib2_overlap halobias_nfw_nopart_${i}_galaxies_18 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_${overlap}_18

elif	[[ $overlap == "d_o" ]]
then
	./IDfib2_double_overlap halobias_nfw_nopart_${i}_galaxies_18 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_${overlap}_18

else
	cp halobias_nfw_nopart_${i}_galaxies_18 halobias_nfw_nopart_${i}_galaxies_${overlap}_18

fi	
	

$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile3 --galaxy_file=halobias_nfw_nopart_${i}_galaxies_${overlap}_18  -output_tag=${output}.$i.$overlap.full_sphere  --theta_max=$maxtheta --n_bins_per_decade=$binning 

done

rm halobias_fof_nfw_${i}_fff_18
#rm halobias_nfw_nopart_${i}_galaxies




done



