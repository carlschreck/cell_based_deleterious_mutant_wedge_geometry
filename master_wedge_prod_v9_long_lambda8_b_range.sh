#!/bin/bash
dir=~/production/deleterious_muts_block_restart_vels_wedge_rotatedBC_initlayer_fast

att=0.0

steps=-1 # don't matter, run until extinct
burnsteps=-1 # don't matter, run until extinct

layerwidth=1.0
propdepth=4.0
bounddepth=4.0

movie=.true.
restart=.true.
bottom=.false.

cutrelax=2000

########################################################
##############     Fig 1 (a) and (b)     ###############
########################################################

Lx=800.0
w0=400.0

s1=-0.06

layerdepth=8.0

for seed in `seq 1 10`
do

b=8e3
dt=2.5e-6
layerskip=400
dataskip=10000
prodskip=50000
restskip=10000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

b=4e3
dt=5e-6
layerskip=800
dataskip=20000
prodskip=100000
restskip=20000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

b=2e3
dt=1e-5
layerskip=1600
dataskip=40000
prodskip=200000
restskip=40000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

b=1e3
dt=2e-5
layerskip=3200
dataskip=80000
prodskip=400000
restskip=80000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

b=5e2
dt=4e-5
layerskip=6400
dataskip=160000
prodskip=800000
restskip=160000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

b=2.5e2
dt=8e-5
layerskip=12800
dataskip=320000
prodskip=1600000
restskip=320000
sbatch $dir/run_cori_v9_long.sh $att $Lx $b $steps $burnsteps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $bottom $s1 $w0 $seed $cutrelax

done
