#!/bin/bash

rbmm=/home/piscioja/makemock/make_random_mock
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
binning=10
maxtheta=0.1
# geometry zspace ifib 
zspace=0
ifib=0
rot=" 0 0 0 "
fgot=" 0.6 " # arbitrary between 0 and 1 for lss_geometry file

#map
map_dir=/data2/jap/SDSSPix/Maps
pix_ns=$map_dir/lss_geometry.dr72.pix
pix_no=$map_dir/lss_geometry_north.dr72.pix
pix_window=$map_dir/window.dr72brvoid0.stripe_trim.pix
pix_overlap=$map_dir/window.dr72brvoid1.overlap.pix
pix_overlap_noweight=$map_dir/window.dr72brvoid1.overlap.noweight.pix
pre="sdssmock_gamma"

sim="consuelo"
lbox=420.0
# cvl_m18 sample : consuelo
# target density : 3.1e-2
samp="main18"

for tag in 10 
do

       npart=`expr $tag \\* 2300000`
       seed=`expr $tag \\* 17`

#	for foot in ns overlap sphere
	for foot in ns 
	do

       if   [[ $foot == "overlap" ]]
       then
               geometry=1
               pix=$pix_overlap_noweight
       elif [[ $foot == "ns" ]]
       then
               geometry=1
               pix=$pix_ns
       else
               geometry=0
               pix="none"
       fi
	
	args=" $geometry $zspace $ifib $rot "  

	$rbmm $lbox $npart $args 0.02 0.042 $pix $fgot $seed  0 > ${pre}_$samp.rand_${tag}x_weighted.rdcz.dat
	
	$stompdirectory/stomp_galaxy_autocorrelation --map_file=$pix --galaxy_file=${pre}_$samp.rand_${tag}x_weighted.rdcz.dat  -output_tag=${pre}_$samp.rand_${tag}x_weighted.$foot.rdcz  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index

	done

done

# cvl_m19 sample : consuelo
# target density : 1.55e-2
samp="main19"

for tag in 10
do

  
        npart=`expr $tag \\* 1150000`
       seed=`expr $tag \\* 19` 

#        for foot in ns overlap sphere
	for foot in ns 
        do

       if   [[ $foot == "overlap" ]]
       then
               geometry=1
               pix=$pix_overlap_noweight
       elif [[ $foot == "ns" ]]
       then
               geometry=1
               pix=$pix_ns
       else
               geometry=0
               pix="none"
       fi

        args=" $geometry $zspace $ifib $rot "

        $rbmm $lbox $npart $args 0.02 0.067 $pix $fgot $seed  0 > ${pre}_$samp.rand_${tag}x_weighted.rdcz.dat

        $stompdirectory/stomp_galaxy_autocorrelation --map_file=$pix --galaxy_file=${pre}_$samp.rand_${tag}x_weighted.rdcz.dat  -output_tag=${pre}_$samp.rand_${tag}x_weighted.$foot.rdcz  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index

        done

done

samp="main20"
sim="esmeralda"
lbox=640.0

for tag in 10
do


        npart=`expr $tag \\* 1630000`
        seed=`expr $tag \\* 20`


#        for foot in ns overlap sphere
	for foot in ns 
        do

       if   [[ $foot == "overlap" ]]
       then
               geometry=1
               pix=$pix_overlap_noweight
       elif [[ $foot == "ns" ]]
       then
               geometry=1
               pix=$pix_ns
       else
               geometry=0
               pix="none"
       fi

        args=" $geometry $zspace $ifib $rot "

        $rbmm $lbox $npart $args 0.02 0.106 $pix $fgot $seed  0 > ${pre}_$samp.rand_${tag}x_weighted.rdcz.dat

        $stompdirectory/stomp_galaxy_autocorrelation --map_file=$pix --galaxy_file=${pre}_$samp.rand_${tag}x_weighted.rdcz.dat  -output_tag=${pre}_$samp.rand_${tag}x_weighted.$foot.rdcz  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index
	done

done

samp="main21"
sim="carmen"
lbox=1000.0

for tag in 10
do


       npart=`expr $tag \\* 1250000`
       seed=`expr $tag \\* 21`



#        for foot in ns overlap sphere
	for foot in ns 
        do

       if   [[ $foot == "overlap" ]]
       then
               geometry=1
               pix=$pix_overlap_noweight
       elif [[ $foot == "ns" ]]
       then
               geometry=1
               pix=$pix_ns
       else
               geometry=0
               pix="none"
       fi

        args=" $geometry $zspace $ifib $rot "

        $rbmm $lbox $npart $args 0.02 0.165 $pix $fgot $seed  0 > ${pre}_$samp.rand_${tag}x_weighted.rdcz.dat
        $stompdirectory/stomp_galaxy_autocorrelation --map_file=$pix --galaxy_file=${pre}_$samp.rand_${tag}x_weighted.rdcz.dat  -output_tag=${pre}_$samp.rand_${tag}x_weighted.$foot.rdcz  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index
        done

done

