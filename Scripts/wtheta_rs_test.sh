#!/bin/bash

#Stomp Parameters
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
maxtheta=.10
binning=10
data_dir=/data2/jap/GalaxyCatalogues/Data
rs_dir=/data2/jap/GalaxyCatalogues/roman_test

#Maps
map_dir=/data2/jap/SDSSPix/Maps
map_sdss=$map_dir/window.dr72brvoid0.stripe_trim.pix
map_lss=$map_dir/lss_geometry.dr72.stripe_trim.pix
map_full=$map_dir/lss_geometry.fullsphere.pix
map_overlap1=$map_dir/window.dr72brvoid1.overlap.balkanized.noweight.pix
map_overlap2=$map_dir/window.dr72brvoid1.overlap.noweight.pix


#for sample in 18 19 20 21
#do
#for fib in fib0 fib1
#do
#$stompdirectory/stomp_galaxy_autocorrelation -map_file=$map_sdss -galaxy_file=vollim_Mr${sample}_$fib.rdcz -output_tag=vollim_Mr${sample}_$fib -n_bins_per_decade=$binning  -theta_max=$maxtheta -single_index
#$stompdirectory/stomp_galaxy_autocorrelation -map_file=$map_overlap2 -galaxy_file=vollim_Mr${sample}_$fib.rdcz -output_tag=vollim_Mr${sample}_$fib.overlap -n_bins_per_decade=$binning  -theta_max=$maxtheta -single_index
#done
#done


#for sample in 18 19 20 21
#do
#for fib in fib0 fib1
#do
#$stompdirectory/stomp_galaxy_autocorrelation_jack -map_file=$map_sdss -galaxy_file=vollim_Mr${sample}_$fib.rdcz -output_tag=vollim_Mr${sample}_$fib.20rand.short -n_bins_per_decade=$binning  -theta_max=$maxtheta -single_index -n_random=20 
#$stompdirectory/stomp_galaxy_autocorrelation_jack -map_file=$map_overlap2 -galaxy_file=vollim_Mr${sample}_$fib.rdcz -output_tag=vollim_Mr${sample}_$fib.20rand.overlap.short -n_bins_per_decade=$binning  -theta_max=$maxtheta -single_index -n_random=20 
#done
#done

for sample in 18 19 20 21
do
/home/piscioja/makemock/trim_geometry 0.6  /data2/jap/SDSSPix/Maps/window.dr72brvoid1.overlap.balkanized.noweight.pix < /data2/jap/GalaxyCatalogues/Data/vollim_Mr${sample}_fib0.rdcz > vollim_Mr${sample}_fib0_overlap.rdcz 
$stompdirectory/stomp_galaxy_autocorrelation_jack -map_file=$map_overlap1 -galaxy_file=vollim_Mr${sample}_fib0_overlap.rdcz -output_tag=vollim_Mr${sample}_fib0_overlap.20rand -n_bins_per_decade=$binning  -theta_max=$maxtheta -single_index -n_random=20

done
