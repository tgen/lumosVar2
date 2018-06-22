#!/bin/bash
#SBATCH -N 1
#SBATCH -J lumosVarMain
#SBATCH -t 14-00:00:0
#SBATCH --mail-type=all
#SBATCH --mail-user=rhalperin@tgen.org
#SBATCH -o lumosVarMain.o%j

module load MCR/9.0
#module load samtools/1.3.1
#export PATH=/home/rhalperin/bin/samtools-1.5/bin:$PATH
#export PERL5LIB=$PERL5LIB:/home/rhalperin/perl5lib/lib/
source /home/rhalperin/multiSampleCaller/gvm/setup.sh

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
