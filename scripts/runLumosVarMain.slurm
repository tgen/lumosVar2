#!/bin/bash
#SBATCH -N 1
#SBATCH -J lumosVarMain
#SBATCH -t 14-00:00:0
#SBATCH --mail-type=all
#SBATCH -o lumosVarMain.o%j

source scripts/setup.sh

TMPDIR=/scratch/rhalperin/TMP_$RANDOM
mkdir $TMPDIR
echo "Running MultiSampleCaller"
echo "CONFIG is ${CONFIG}"
export MCR_CACHE_ROOT=$TMPDIR
export TMP=$TMPDIR
echo "Cache going to $TMPDIR"

date
IFS='. ' read -a NAMES <<< $CONFIG
./bin/lumosVarMain ${CONFIG} $SLURM_CPUS_ON_NODE >${NAMES[0]}_LumosVarOut.txt
date

rm -rf $TMPDIR
