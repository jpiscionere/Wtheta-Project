#! /bin/bash

output=stomp_sdssmock_main20_esmeralda
maxtheta=0.1
binning=10
#fof=/home/piscioja/halobias/halobias_fof_nfw
fof=halobias_fof_nfw_quiet
so=/home/piscioja/halobias/halobias_so_nfw2
map_dir=/home/piscioja/SDSSPix/Maps
map_sdss=$map_dir/window.dr72brvoid0.stripe_trim.pix
map_lss=$map_dir/lss_geometry.dr72.stripe_trim.pix
map_full=$map_dir/lss_geometry.fullsphere.pix
map_overlap1=$map_dir/window.dr72brvoid1.overlap.balkanized.noweight.pix
map_overlap2=$map_dir/window.dr72brvoid1.overlap.noweight.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=11.8067544233
siglogM=0.0963472124455
logM0=12.6369512043
logM1=13.1533827128
alpha=1.06123071572
delta_vir=377 #linking length of 0.2
Mstar=2.29E12
redshift=0.082

logMmin=11.5474784286 
siglogM=0.0980446986406 
logM0=11.1088487112 
logM1=13.0620572299 
alpha=1.16251908367 
gamma=0.241600544769 
fgal=1.8570672102

logMmin=11.7565254961 
siglogM=0.0972729443771 
logM0=12.0985139989 
logM1=13.1786228725 
alpha=1.08475003121 
gamma=0.205950610178 
fgal=1.81881523299 


logMmin=11.7577328298   
siglogM=0.0972008997931 
logM0=12.7430056861   
logM1=13.0177787909   
alpha=1.04953050812   
gamma=0.204769578295  
fgal=1.7682103135   


for i in $(seq 3011 1 3012) 
do


#for gamma in 0.206254091967 
#do

#for fgal in 1.79335961349 
#do



#time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < /ssd1/Research/halo_files/Esmeralda/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > halobias_fof_nfw_${i}_fff

time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < /ssd1/Research/halo_files/Esmeralda/Fof/Esmeralda_${i}_z0p000_fof_b0p2.dpp.halos > halobias_fof_nfw_${i}_fff

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.106 $map_lss 0.6 0 < halobias_fof_nfw_${i}_fff > halobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation_full --map_file=$map_lss --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.nfw.fof.lss_geometry2  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 

#time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$map_overlap1 --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.nfw.fof.overlap2  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index

rm halobias*


done
#done
#done



