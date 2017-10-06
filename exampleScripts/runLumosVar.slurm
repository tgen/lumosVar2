#!/usr/bin/env bash

#SBATCH --job-name=lumosVar
#SBATCH --time=48:00:00
#SBATCH --output=lumosVar.o%j
#SBATCH --nodes=1

CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)


cd $SLURM_SUBMIT_DIR

#LIBRARY_PATH=/packages/gcc/5.1.0/lib64:/packages/gcc/5.1.0/lib:/home/skulkarni/packages/htslib/git-head/lib:/home/skulkarni/packages/yaml/0.1.7/lib:/cm/shared/apps/slurm/17.02.3/lib64/slurm:/cm/shared/apps/slurm/17.02.3/lib64
#LD_LIBRARY_PATH=/packages/gcc/5.1.0/lib64:/packages/gcc/5.1.0/lib:/home/skulkarni/packages/htslib/git-head/lib:/home/skulkarni/packages/yaml/0.1.7/lib:/cm/shared/apps/slurm/17.02.3/lib64/slurm:/cm/shared/apps/slurm/17.02.3/lib64
#PATH=/packages/gcc/5.1.0/bin:/home/skulkarni/packages/htslib/git-head/bin:/cm/shared/apps/slurm/17.02.3/sbin:/cm/shared/apps/slurm/17.02.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/sbin:/cm/local/apps/environment-modules/3.2.10/bin:/opt/dell/srvadmin/bin
#C_INCLUDE_PATH=/packages/gcc/5.1.0/include:/home/skulkarni/packages/htslib/git-head/include:/home/skulkarni/packages/yaml/0.1.7/include
#_LMFILES_=/packages/modules-3.2.9/modulefiles/slurm/17.02.3:/home/skulkarni/privatemodules/yaml/0.1.7:/home/skulkarni/privatemodules/htslib/git-head:/packages/modules-3.2.9/modulefiles/gcc/5.1.0:/packages/modules-3.2.9/modulefiles//slurm/17.02.3

export MODULEPATH=$MODULEPATH:/home/skulkarni/privatemodules/
source /home/rhalperin/multiSampleCaller/gvm/setup.sh
module load yaml
module load samtools
module load htslib
module load MCR/9.0

GVM_EXEC=/home/rhalperin/multiSampleCaller/gvm/src/gvm

TMPDIR=/scratch/rhalperin/TMP_$RANDOM
mkdir $TMPDIR
date
echo "Running MultiSampleCaller"
echo "CONFIG is ${CONFIG}"
export MCR_CACHE_ROOT=$TMPDIR
export TMP=$TMPDIR
echo "Cache going to $TMPDIR"
echo "GVM: $GVM_EXEC"

#env >slurm_env.txt

#srun -n1 --export=ALL $GVM_EXEC $CONFIG 22

IFS='. ' read -a NAMES <<< $CONFIG
FAIL_COUNT=0
pids=""
echo "### Submitting to queue to run lumosVar on ${NAMES[0]}"
for CHR in "${CHRS[@]}"
do
        echo "submiting step ${CHR}"
        #touch ${NAMES[0]}_Step${CHR}.lumosVarInQueue
        srun -n1 --export=ALL $GVM_EXEC --conf $CONFIG --chr $CHR &
        pids+=" $!"
        #if [ $? -eq 0 ] ; then
        #	touch ${NAMES[0]}_Step${CHR}.lumosVarPass
	#else
        #	touch ${NAMES[0]}_Step${CHR}.lumosVarFail
	#	(( FAIL_COUNT++ ))
	#fi
	#rm ${NAMES[0]}_Step${CHR}.lumosVarInQueue
done
#wait
for p in $pids; do
        if wait $p; then
                echo "Process $p success"
        else
                echo "Process $p fail"
		(( FAIL_COUNT++ ))
        fi
done

echo "### Finished to gvm with $FAIL_COUNT fails"
date

if [ $FAIL_COUNT -eq 0 ] ; then
	touch ${NAMES[0]}_Main.lumosVarInQueue
        echo "### Starting lumosVarMain on $SLURM_CPUS_ON_NODE nodes"
	srun --export=ALL -n 1 -c $SLURM_CPUS_ON_NODE  bin/lumosVarMain ${CONFIG} >${NAMES[0]}_LumosVarOut.txt
        if [ $? -eq 0 ] ; then
                touch ${NAMES[0]}_Main.lumosVarPass
                exit 0
        else
               echo "### lumosVarMain failed" 
               touch ${NAMES[0]}_Main.lumosVarFail
               exit 1
        fi
	date
else
	exit 1
fi

rm -rf $TMPDIR
