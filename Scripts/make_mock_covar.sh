#!/bin/bash

source_dir=/home/piscioja/Clustering/WthetaPaper/Source
data_dir=/home/piscioja/Clustering/WthetaPaper/Data/Wtheta
mock_dir=/data2/jap/GalaxyCatalogues/NFW
covar_dir=/home/piscioja/Clustering/WthetaPaper/Data/Mock_Covar
rand_dir=/data2/jap/GalaxyCatalogues/Randoms
cov2=$source_dir/covariance_fitting2

bins2=20 #number of total bins
print_var=0 #print the variance instead of the covariance matrix
diag_error=1 #use diagonal errors
slope0=-1.4
d_s=0.00075
n_s=2000
amp0=0.0
d_amp=0.01
n_amp=2000

args="$bins2 $print_var $diag_error $slope0 $d_s $n_s $amp0 $d_amp $n_amp"


for mr in 20 21
do
	pre=vollim_Mr${mr}

	if [[ $print_var == 1 ]]
	then
		post=variance
	else
		post=mock_covar
	fi

	for overlap in 0 1
	do

		if [[ $overlap == 0 ]]
		then 
			data=$data_dir/Wtheta_vollim_Mr${mr}_fib0.20rand.short
			rand=$rand_dir/Wtheta_sdssmock_gamma_main${mr}.rand_10x_weighted.ns.rdcz.wtheta
			models=$mock_dir/$mr/Wtheta*lss_geometry.wtheta
			tag3=fib0
		else
			data=$data_dir/Wtheta_vollim_Mr${mr}_fib0.20rand.overlap.short
			rand=$rand_dir/Wtheta_sdssmock_gamma_main${mr}.rand_10x_weighted.overlap.rdcz.wtheta
			models=$mock_dir/$mr/Wtheta*overlap.wtheta
			tag3=fib0_overlap
		fi
	
		awk '{if(NR <= '$bins2'){print $1,$2,$3}}'<$data >data_in.tmp	
		awk '{if(NR <= '$bins2'){print $3}}'<$rand >rand_in.tmp	
		$cov2 data_in.tmp rand_in.tmp $args $models > Covariance_${pre}_${tag3}.$post 

	done
done
