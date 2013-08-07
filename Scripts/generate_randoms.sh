#!/bin/bash

rb=/home/mcbridck/code/lasdamas/bin/randbox_fastfood
mm=/home/mcbridck/code/lasdamas/bin/makemock
rbmm=/home/piscioja/makemock/make_random_mock
astrostomp=/home/piscioja/astrostomp/astro-stomp-read-only/examples/stomp_galaxy_autocorrelation
mapfile=/home/piscioja/GalaxyCatalogues/window.dr72brvoid0.stripe_trim.pix
data=/home/piscioja/Clustering/Data
thetamax=0.01
n_jack=20
nbins=10
# geometry zspace ifib 
zspace=0
ifib=0
rot=" 0 0 0 "
fgot=" 0.6 " # arbitrary between 0 and 1 for lss_geometry file
# full sphere

outdir="randoms"
pre="sdssmock_gamma"



pix_ns="/data0/vagc-nyu/sdsspix/dr72/lss_geometry.dr72.pix"
pix_no="/data0/vagc-nyu/sdsspix/dr72/lss_geometry_north.dr72.pix"
pix_window="/data2/jap/GalaxyCatalogues/window.dr72brvoid0.stripe_trim.pix"
pix_overlap="/data2/jap/SDSSPix/overlap_pixel_balkanized_short.noweight.pix"
sim="consuelo"
lbox=420.0
# cvl_m18 sample : consuelo
# target density : 3.1e-2
#samp="main18"


#for tag in 10 20 30 40
#for tag in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 

#do
	
#	npart=`expr $tag \\* 2300000`
#	npart=46000000
#	seed=`expr $tag \\* 17`

#	$rb $lbox $npart $seed > $samp.$sim.rand_${tag}x.ff 
#	$rb $lbox $npart $seed > $samp.$sim.rand_20x.ff 

#	for foot in ns no sphere
#	for foot in no
#
#	do
#   	 
#	if   [[ $foot == "no" ]]
#   	then
#        	geometry=1
#        	pix=$pix_no
#    	elif [[ $foot == "ns" ]]
#    	then
#        	geometry=1
#        	pix=$pix_ns
#    	else
#        	geometry=0
#        	pix="none"
#    	fi
#    	args=" $geometry $zspace $ifib $rot "	
#
#		$mm $args 0.02 0.042 $pix $fgot 0 < $samp.$sim.rand_${tag}x.ff > ${pre}_$samp.rand_${tag}x_${foot}.rdcz.dat	
#		$mm $args 0.02 0.042 $pix $fgot 0 < $samp.$sim.rand_20x.ff > ${pre}_$samp.rand_20x_${foot}_${tag}.rdcz.dat	
#
#	done
#
#	rm $samp.$sim.rand_20x.ff	
#	
#	$astrostomp --map_file=$mapfile --galaxy_file=$data/vollim_Mr18_fib0.rdcz --randoms_file=${pre}_${samp}.rand_20x_no_${tag}.rdcz.dat -single_index -output_tag=${samp}_fib0_20x_randoms_${tag} --n_bins_per_decade=$nbins --theta_max=$thetamax --n_jackknife=$n_jack
#	$astrostomp --map_file=$mapfile --galaxy_file=$data/vollim_Mr18_fib0.rdcz --randoms_file=${pre}_${samp}.rand_20x_no_${tag}.rdcz.dat -single_index -output_tag=${samp}_fib0_20x_randoms_${tag} --n_bins_per_decade=$nbins --theta_max=$thetamax 
#	rm ${pre}_$samp.rand_20x_${foot}_${tag}.rdcz.dat
#done




# cvl_m19 sample : consuelo
# target density : 1.55e-2
samp="main19"

#for tag in 10 20 30 40 
#for tag in 20
#do

#        npart=`expr $tag \\* 1150000`
#	seed=`expr $tag \\* 19`

#        $rb $lbox $npart $seed > $samp.$sim.rand_${tag}x.ff        

#        for foot in ns no sphere
#        for foot in sphere 
#        do
#         
#        if   [[ $foot == "no" ]] 
#        then
#                geometry=1
#                pix=$pix_no
#        elif [[ $foot == "ns" ]]
#        then
#                geometry=1
#                pix=$pix_ns
#        else
#                geometry=0
#                pix="none"
#        fi
#        args=" $geometry $zspace $ifib $rot "   
#        
#                $mm $args 0.02 0.067 $pix $fgot 0 < $samp.$sim.rand_${tag}x.ff > ${pre}_$samp.rand_${tag}x_${foot}.rdcz.dat                
#       
#        done
        
#	rm $samp.$sim.rand_${tag}x.ff
#	$astrostomp --map_file=$mapfile --galaxy_file=$data/vollim_Mr19_fib0.rdcz --randoms_file=${pre}_${samp}.rand_${tag}x_no.rdcz.dat -single_index -output_tag=${samp}_fib0_${tag}x_randoms --n_bins_per_decade=$nbins --theta_max=$thetamax --n_jackknife=$n_jack
#
#done
#
#

samp="main20"
sim="esmeralda"
lbox=640.0


#for tag in 10 20 30 40 
for tag in 10
do

        npart=`expr $tag \\* 1630000`
	seed=`expr $tag \\* 20`

#        $rb $lbox $npart $seed > $samp.$sim.rand_${tag}x.ff

#        for foot in ns no sphere
        for foot in ns 
        do

        if   [[ $foot == "no" ]]
        then
                geometry=1
                pix=$pix_no
        elif [[ $foot == "ns" ]]
        then
                geometry=1
                pix=$pix_ns
        else
                geometry=0
                pix="none"
        fi
        args=" $geometry $zspace $ifib $rot "

#                $mm $args 0.02 0.106 $pix_window 0.6 0 < $samp.$sim.rand_${tag}x.ff > ${pre}_$samp.rand_${tag}x_ns.rdcz.dat
#                $mm $args 0.02 0.106 $pix_overlap 0.00001 0 < $samp.$sim.rand_${tag}x.ff > ${pre}_$samp.rand_${tag}x_overlap.rdcz.dat
                $rbmm $lbox $npart $args 0.02 0.106 $pix_window $fgot $seed  0 > ${pre}_$samp.rand_${tag}x_weighted.rdcz.dat
     #           $rbmm $lbox $npart $args 0.02 0.106 $pix_ns $fgot $seed  0 > sdss_geo_randoms 
                #valgrind -v --leak-check=full --show-reachable=yes --track-origins=yes $rbmm $lbox 10 $args 0.02 0.106 $pix $fgot $seed  0 > test 
                #valgrind -v  $rbmm $lbox 10 $args 0.02 0.106 $pix $fgot $seed  0 > test 
       
        done
        
#	rm $samp.$sim.rand_${tag}x.ff	
#	$astrostomp --map_file=$mapfile --galaxy_file=$data/vollim_Mr20_fib0.rdcz --randoms_file=${pre}_${samp}.rand_${tag}x_no.rdcz.dat -single_index -output_tag=${samp}_fib0_${tag}x_randoms --n_bins_per_decade=$nbins --theta_max=$thetamax --n_jackknife=$n_jack
#
done



sim="carmen"
lbox=1000.0
samp="main21"


#for tag in 10 20 30 40 
#for tag in 10 
#do

#        npart=`expr $tag \\* 1250000`
#	seed=`expr $tag \\* 21`
#        $rb $lbox $npart $seed > $samp.$sim.rand_${tag}x.ff

#         for foot in ns no sphere
#        for foot in ns 
#        do

#        if   [[ $foot == "no" ]]
#        then
#                geometry=1
#                pix=$pix_no
#        elif [[ $foot == "ns" ]]
#        then
#                geometry=1
#                pix=$pix_ns
#        else
#                geometry=0
#                pix="none"
#        fi
#        args=" $geometry $zspace $ifib $rot "
#
#                $mm $args 0.02 0.165 $pix $fgot 0 < $samp.$sim.rand_${tag}x.ff > ${pre}_$samp.rand_${tag}x_${foot}.rdcz.dat
#       
#        done
#        
#	rm $samp.$sim.rand_${tag}x.ff
#	$astrostomp --map_file=$mapfile --galaxy_file=$data/vollim_Mr21_fib0.rdcz --randoms_file=${pre}_${samp}.rand_${tag}x_no.rdcz.dat -single_index -output_tag=${samp}_fib0_${tag}x_randoms --n_bins_per_decade=$nbins --theta_max=$thetamax

#done



