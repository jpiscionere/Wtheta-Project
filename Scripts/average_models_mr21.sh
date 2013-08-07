#! /bin/bash

sim=Carmen
model_dir=/data2/jap/GalaxyCatalogues/FgalGamma/$sim
output=Wtheta_main21_carmen
input=Wtheta_stomp_sdssmock_main21_carmen

s_box=2020
n_boxes=4



for gamma1 in $(seq 1.5 0.02 2.5 )
do

for fgal1 in $(seq 0.01 0.01 0.51)
do

gamma="$(printf "%4.3f" "$gamma1")"
fgal="$(printf "%4.3f" "$fgal1")"


#Clearing total.dat
awk '{print $1/$1 - 1}' <$model_dir/total.dat>$model_dir/tmp1.dat
mv $model_dir/tmp1.dat $model_dir/total.dat

for (( i = $s_box; i < $s_box + $n_boxes ; i++))
do


data_file=$model_dir/$i/${input}.$i.f$fgal.g${gamma}.fof.lss_geometry.wtheta

awk '{DD = $3; printf("%6.5e\n",DD)}' <$data_file >$model_dir/tmp.dat
pr -m -t -s $model_dir/total.dat $model_dir/tmp.dat | awk '{sum=$1+$2; printf("%6.5e\n",sum)}' >$model_dir/tmp2.dat
mv $model_dir/tmp2.dat $model_dir/total.dat
rm $model_dir/tmp.dat

done


awk '{average = $1/'$n_boxes'; printf("%6.5e\n",average)}' < $model_dir/total.dat > $model_dir/Averages/${output}.average.f$fgal.g${gamma}.wtheta
echo " '$fgal' '$gamma'"



done
done
