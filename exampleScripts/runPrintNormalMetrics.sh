#!/bin/bash -l
CONFIG=$1


STEP=0
STEP_COUNT=24
#STEP_COUNT=`ls $chrList/*list | wc -l`
IFS='. ' read -a NAMES <<< $CONFIG
echo "### Submitting to queue to run lumosVar on ${NAMES[0]}"
while [ ${STEP} -lt $STEP_COUNT ]
do
        (( STEP++ ))
        if [[ -e ${NAMES[0]}_Step${STEP}.lumosVarInQueue || -e ${NAMES[0]}_Step${STEP}.lumosVarPass || -e ${NAMES[0]}_Step${STEP}.lumosVarFail ]] ; then
        	echo "### LumosVar is already done, failed, or inqueue for ${CONFIG} step ${STEP}"
                continue
        fi
	echo "submitting step ${STEP}"
	touch ${NAMES[0]}_Step${STEP}.lumosVarInQueue
	qsub -v CONFIG=${CONFIG},CHR=${STEP} /home/rhalperin/pbs_scripts/printNormalMetricsByChr.pbs
	STEP_DONE[$STEP]=0
	sleep 2
done

        	
