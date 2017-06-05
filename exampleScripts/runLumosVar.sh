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
	echo "submiting step ${STEP}"
	qsub -v CONFIG=${CONFIG},CHR=${STEP} /home/rhalperin/pbs_scripts/runLumosVarPreprocess.pbs
	touch ${NAMES[0]}_Step${STEP}.lumosVarInQueue
	STEP_DONE[$STEP]=0
	sleep 2
done

PASS_COUNT=0
until [ $PASS_COUNT -eq $STEP_COUNT ]
do
	STEP=0
	while [ ${STEP} -lt $STEP_COUNT ]
 	do
		(( STEP++ ))
		if [[ ${STEP_DONE[$STEP]} -ne 1 && -e ${NAMES[0]}_Step${STEP}.lumosVarPass ]]; then
			STEP_DONE[$STEP]=1
			(( PASS_COUNT++ ))
		fi				
	done
	if [[ $PASS_COUNT -lt $STEP_COUNT ]]; then
		echo "current pass count ${PASS_COUNT}"
		sleep 300
	fi 
done

echo "running main LumosVar"
qsub -v CONFIG=${CONFIG} /home/rhalperin/pbs_scripts/runLumosVarMain.pbs
        	
