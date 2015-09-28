#!/bin/sh 

for x in 4e-2
do
for vy in 0
do
for ene in 4.823e5
do
for sol_length in 0.3
do
for nstep in 6001
do

rm *dat

sed s/INIT_AMP/"$x"/g  sol_scipy_sample.py > sol_scipyc1.py

sed s/INIT_VELO/"$vy"/g  sol_scipyc1.py > sol_scipyc2.py

sed s/KINE_ENE/"$ene"/g  sol_scipyc2.py > sol_scipyc3.py

sed s/SOL_LENGTH/"$sol_length"/g  sol_scipyc3.py > sol_scipyc4.py

sed s/NSTEP/"$nstep"/g  sol_scipyc4.py > sol_scipy.py

/usr/local/bin/python sol_scipy.py

rm sol_scipyc*py

mv alldata.dat alldata_amp"$x"_velo"$vy"_kene"$ene"_sleng"$sol_length"_step"$nstep".dat

done
done
done
done
done
