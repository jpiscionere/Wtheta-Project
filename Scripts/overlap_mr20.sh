#! /bin/bash
fof=/home/piscioja/halobias/halobias_fof_nfw
output=stomp_sdssmock_main20_esmeralda
map_dir=/data2/jap/SDSSPix/Maps
maxtheta=10
binning=10
mapfile1=$map_dir/lss_geometry.dr72.stripe_trim.pix
mapfile3=$map_dir/lss_geometry.fullsphere.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
gamma=2.0
fgal=0.2
Mstar=2.29E12
redshift=0.082

#$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile1 --galaxy_file=/data2/jap/GalaxyCatalogues/Randoms/sdssmock_gamma_main21.rand_10x_ns.rdcz.dat  -output_tag=sdssmock_gamma_main21.RR  --theta_max=$maxtheta --n_bins_per_decade=$binning


for i in $(seq 3001 1 3010)
do



#/home/piscioja/halobias/halobias_nfw_nopart_fgal 3 4 1 11.8013669053163 0.0973817154668221 12.687852882335 13.0400932413139 1.05768884410071 1 1 $gamma $fgal trash  -$i /data0/LasDamas/Esmeralda/$i/Esmeralda_${i}_z0p082_fof_b0p2.000*.bgc < /data0/LasDamas/Esmeralda/$i/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > halobias_nfw_nopart_${i}_fff

time $fof 3 4 1 11.8013669053163 0.0973817154668221 12.687852882335 13.0400932413139 1.05768884410071 1 $gamma $fgal -20.0 200.0 $Mstar $redshift trash -$i < /data0/LasDamas/Esmeralda/$i/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > halobias_nfw_nopart_${i}_fff


/home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.106 /data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix 0.6 0 < halobias_nfw_nopart_${i}_fff > halobias_nfw_nopart_${i}_galaxies

for overlap in n_o s_o d_o n_c 


do

if	[[ $overlap == "n_o" ]]
then
	./IDfib2 halobias_nfw_nopart_${i}_galaxies 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_$overlap

elif	[[ $overlap == "s_o" ]]
then
	./IDfib2_overlap halobias_nfw_nopart_${i}_galaxies 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_${overlap}

elif	[[ $overlap == "d_o" ]]
then
	./IDfib2_double_overlap halobias_nfw_nopart_${i}_galaxies 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_$overlap

elif	[[ $overlap == "j_p" ]]
then
	./IDfib2_justpairs halobias_nfw_nopart_${i}_galaxies 1 -$i 0 > halobias_nfw_nopart_${i}_galaxies_$overlap
	
else
	cp halobias_nfw_nopart_${i}_galaxies halobias_nfw_nopart_${i}_galaxies_$overlap

fi	
	

$stompdirectory/stomp_galaxy_autocorrelation --map_file=$mapfile3 --galaxy_file=halobias_nfw_nopart_${i}_galaxies_${overlap}  -output_tag=${output}.$i.$overlap.fgal_gamma  --theta_max=$maxtheta --n_bins_per_decade=$binning

done

rm halobias_nfw_nopart_${i}_fff
#rm halobias_nfw_nopart_${i}_galaxies




done



