#!/usr/bin/env bash

#SBATCH --job-name=normalMetrics
#SBATCH --output=normalMetrics.o%j
#SBATCH --array=0-23
#SBATCH --time=3-00:00:00

CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
echo "task is $SLURM_ARRAY_TASK_ID"

source scripts/setup.sh
echo "config is $CONFIG"
echo "chr is $CHR"

python scripts/runNormalMetrics.py $CONFIG $CHR


