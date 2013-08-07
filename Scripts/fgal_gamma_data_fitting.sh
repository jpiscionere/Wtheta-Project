use_mocks=0

s_f=0 #covariance matrix input format
i_f=1 #invert the matrix?
diag_err=0 #use diagonal errors
u_v=0 #use variance instead of jackknife errors?


bins2=20 #number of total bins
s_bin=0 #first bin you care about
e_bin=20 #last bin you care about
g0=1.5 #gamma0
d_g=0.008 #delta_gamma
n_g=51 #number of gamma steps
f0=0.010 #fgal0
d_f=0.005 #delta_fgal
n_f=41 #number of fgal steps

args="$s_f $i_f $u_v $diag_err $bins2 $s_bin $e_bin $g0 $d_g $n_g $f0 $d_f $n_f"

main=/home/piscioja/Clustering/WthetaPaper
data_dir=$main/Data/Wtheta
out_dir=/hd0/Research/Clustering/FgalGammaOutputs
source_dir=$main/Source
rand_dir=/hd0/Research/Clustering/Randoms
cov1=covariance_fitting
cov2=covariance_fitting_sigma
cov3=covariance_fitting_fgalgamma
model_dir=/home/piscioja/Clustering/WthetaPaper/Data/Average_Models

for mr in 20
do 
	if [[ $mr == 20 ]]
	then
		m_dir=Esmeralda
		tag2=main20_esmeralda
	else
		m_dir=Carmen
		tag2=main21_carmen
	fi
	
	for overlap in 0 1
	do

		if [[ $overlap == 0 ]]
		then
			data=Wtheta_vollim_Mr${mr}_fib0.20rand.short
			rand=$rand_dir/Wtheta_sdssmock_gamma_main${mr}.rand_10x_weighted.ns.rdcz.wtheta
			
			if [[ $use_mocks == 1 ]]
                        then
                                variance=Covariance_vollim_Mr${mr}_fib0.variance
                                covar=Covariance_vollim_Mr${mr}_fib0.mock_covar
                                tag4=mock_covar
                        else
                                covar=$data_dir/CovarWtheta_vollim_Mr${mr}_fib0.20rand.short
                                variance=variance_in.tmp
                                tag4=jack_covar
                        fi
	
			
			
			tag3=fib0	
		else

			data=Wtheta_vollim_Mr${mr}_fib0.20rand.overlap.short
			rand=$rand_dir/Wtheta_sdssmock_gamma_main${mr}.rand_10x_weighted.overlap.rdcz.wtheta
	
			if [[ $use_mocks == 1 ]]
                        then
                                variance=Covariance_vollim_Mr${mr}_fib0_overlap.variance
                                covar=Covariance_vollim_Mr${mr}_fib0_overlap.mock_covar
                                tag4=mock_covar
                        else
                                covar=$data_dir/CovarWtheta_vollim_Mr${mr}_fib0.20rand.overlap.short
                                variance=variance_in.tmp
                                tag4=jack_covar
                        fi
		

			tag3=fib0_overlap
		fi
	
		
		awk '{if((NR >= '$s_bin') && (NR <= '$e_bin')){print $1,$2,$3}}'<$data_dir/$data >data_in.tmp
		awk '{if((NR >= '$s_bin') && (NR <= '$e_bin')){print $1}}'<$variance >variance_in.tmp
		awk '{if((NR >= '$s_bin') && (NR <= '$e_bin')){print $3}}'<$rand >rand_in.tmp


		$source_dir/$cov3 data_in.tmp variance_in.tmp rand_in.tmp $args $model_dir/$m_dir  $tag2  < $covar > $out_dir/Covariance_fgal_gamma_${tag2}_${tag3}.$tag4.chisquared


	done
done



