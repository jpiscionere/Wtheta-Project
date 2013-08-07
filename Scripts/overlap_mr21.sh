#! /bin/bash
fof=/home/piscioja/halobias/halobias_fof_nfw
output=stomp_sdssmock_main21_carmen
maxtheta=10
binning=10
mapfile2=/data0/cameron/for_jen/Archive/lss_geometry.dr72.stripe_trim.pix
mapfile3=/data2/jap/GalaxyCatalogues/lss_geometry.fullsphere.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
gamma=1
fgal=1
Mstar=1.97E12
redshift=0.132

for i in $(seq 2001 1 2010)
do

#/home/piscioja/halobias/halobias_nfw_nopart_fgal 3 4 1 12.76422334 0.71194502 13.07148775 13.70323156 1.07455104 1 1 $gamma $fgal trash  -$i /data0/LasDamas/Carmen/$i/Carmen_${i}_z0p132_fof_b0p2.00**.bgc < /data0/LasDamas/Carmen/$i/Carmen_${i}_z0p132_fof_b0p2.dpp.halos > halobias_nfw_${i}_fff

#time $fof 3 4 1 12.76422334 0.71194502 13.07148775 13.70323156 1.07455104 1 $gamma $fgal -21.0 200.0 $Mstar $redshift trash -$i < /data0/LasDamas/Carmen/$i/Carmen_${i}_z0p132_fof_b0p2.dpp.halos > halobias_nfw_${i}_fff

#/home/piscioja/halobias/makemock 0 1 0 0 0 0 0.02 0.165 /data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix 0.6 0 < halobias_nfw_${i}_fff > halobias_nfw_${i}_galaxies


for overlap in n_o s_o d_o n_c

do

if      [[ $overlap == "n_o" ]]
then
        ./IDfib2 halobias_nfw_${i}_galaxies 1 -$i 0 > halobias_nfw_${i}_galaxies_$overlap

elif    [[ $overlap == "s_o" ]]
then
        ./IDfib2_overlap halobias_nfw_${i}_galaxies 1 -$i 0 > halobias_nfw_${i}_galaxies_$overlap

elif    [[ $overlap == "d_o" ]]
then
        ./IDfib2_double_overlap halobias_nfw_${i}_galaxies 1 -$i 0 > halobias_nfw_${i}_galaxies_$overlap
else
	cp halobias_nfw_${i}_galaxies halobias_nfw_${i}_galaxies_$overlap

fi


$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile3 --galaxy_file=halobias_nfw_${i}_galaxies_$overlap  -output_tag=$output.$i.$overlap.full_sky  --theta_max=$maxtheta --n_bins_per_decade=$binning 
#$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile2 --galaxy_file=halobias_nfw_${i}_galaxies_$overlap  -output_tag=$output.$i.$overlap  --theta_max=$maxtheta --n_bins_per_decade=$binning 


#rm halobias_nfw_${i}_fff
#rm halobias_nfw_${i}_galaxies

done

done

