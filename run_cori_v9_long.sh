#!/bin/bash

#SBATCH -p shared
#SBATCH -C haswell
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J LRE_VELS

# Set directories
rundir=~/production/deleterious_muts_block_restart_vels_wedge_rotatedBC_initlayer_fast
outdir=/global/cscratch1/sd/cschreck/deleterious_muts_block_restart_vels_wedge_rotatedBC_initlayer_fast

# geometric parameters
ar1=1.01
ar2=2.0
Lx=$2
att=$1

# rates
rate0=1d0
b=$3

# steps
steps=$4
burnsteps=$5
layerskip=$6
dataskip=$7
prodskip=$8
restskip=$9
dt=${10}

# growth layer widths
layerwidth=${11}
layerdepth=${12}
propdepth=${13}
bounddepth=${14}

# run parameters
desync=0.4
seed=-${20}

# logical variables
movie=${15}
restart=${16}
bottom=${17}

# rescue variables
s1=${18}
w0=${19}

# cut steps
cutrelax=${21}

# output files
suffix=v9_L${Lx}_layer${layerdepth}_desync${desync}_b${b}_att${att}_s${s1}_winit${w0}_seed${seed}.dat
prodfile=prod_wedge_$suffix
restfile=restart_wedge_$suffix
enerfile=V_wedge_$suffix
contfile=contacts_wedge_$suffix

cd $outdir

echo $cutsteps1
echo $cutsteps2

# run program
time $rundir/dimers_damped_linear_front_removecells_block_vels_wedge_rotatedBC_initlayer_v9_extinction.o <<EOF
  $ar1
  $ar2
  $Lx
  $att
  $rate0
  $b
  $steps
  $burnsteps
  $layerskip
  $dataskip
  $prodskip
  $restskip
  $dt
  $layerwidth
  $layerdepth
  $propdepth
  $bounddepth
  $desync
  $seed
  $prodfile 
  $restfile 
  $contfile 
  $enerfile 
  $movie
  $restart
  $bottom 
  $s1
  $w0
  $cutrelax
EOF
