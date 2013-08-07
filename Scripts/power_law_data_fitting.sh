#!/bin/bash

use_mocks=0
use_var=0 #use variance instead of jackknife errors?
s_f=1 #covariance matrix input format
i_f=1 #invert the matrix?
diag_error=0 #use diagonal errors

bins2=20 #number of total bins
s_bin=3 #first bin you care about
e_bin=17 #last bin you care about
slope0=-1.4
d_s=0.00075
n_s=2000
amp0=0.0
d_amp=0.01
n_amp=2000

args1="$s_f $i_f $use_var $diag_error $bins2"
args2="$slope0 $d_s $n_s $amp0 $d_amp $n_amp"

main=/home/piscioja/Clustering/WthetaPaper
data_dir=$main/Data/Wtheta
out_dir=$main/Output
source_dir=$main/Source
cov1=$source_dir/covariance_fitting
cov2=$source_dir/covariance_fitting_sigma
cov3=$source_dir/covariance_fitting_fgalgamma


for mr in 20  
do
	for overlap in 0
	do
		pre=vollim_Mr${mr}

                if [[ $overlap == 0 ]]
                then
                        data=$data_dir/Wtheta_vollim_Mr${mr}_fib0.20rand.short
                        out_tag=vollim_Mr${mr}_fib0.20rand

                        if [[ $use_mocks == 1 ]]
			then
				variance=Covariance_vollim_Mr${mr}_fib0.variance
				covar=Covariance_vollim_Mr${mr}_fib0.mock_covar
                        	tag2=mock_covar
			else
                        	covar=$data_dir/CovarWtheta_vollim_Mr${mr}_fib0.20rand.short
				variance=variance_in.tmp
				tag2=jack_covar
			fi

			tag3=fib0
                else

                        data=$data_dir/Wtheta_vollim_Mr${mr}_fib0.20rand.overlap.short
                        out_tag=vollim_Mr${mr}_fib0.20rand.overlap

                        if [[ $use_mocks == 1 ]]
			then
				variance=Covariance_vollim_Mr${mr}_fib0_overlap.variance
                                covar=Covariance_vollim_Mr${mr}_fib0_overlap.mock_covar
				tag2=mock_covar	
        		else
       		        	covar=$data_dir/CovarWtheta_vollim_Mr${mr}_fib0O.20rand.overlap.short
				variance=variance_in.tmp
				tag2=jack_covar              
			fi

			tag3=fib0_overlap
                fi
	


		echo "${pre}_${tag3}"	
		echo "$covar"
		           	     
			awk '{if((NR > '$s_bin') && (NR <= '$e_bin')){print $1,$2,$3}}'<$data >data_tmp
                	awk '{if((NR > '$s_bin') && (NR <= '$e_bin')){print $1}}'<$variance >variance_in.tmp
	

			
			$cov1 data_tmp variance_in.tmp $args1 $s_bin $e_bin $args2 < $covar > $out_dir/Covariance.power_law.${pre}_${tag3}.chisquared	
#			$cov2 $n_amp $n_s $amp0 $d_amp $slope0 $d_s 1 < $out_dir/Covariance.power_law.${pre}_${tag3}_var.chisquared >$out_dir/Covariance.power_law.${pre}_${tag3}.pdf 
#			awk '{if(NR==1) {print "'$mr'","'$overlap'","mock",$1,$2,$3}}'<$out_dir/Covariance.power_law.${pre}_${tag3}.pdf >>outfile


	done

done
