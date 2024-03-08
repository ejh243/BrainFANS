#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/learnModel-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/learnModel-%A_%a.e
#SBATCH --job-name=learnModel

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u


## needs to be executed from the scripts folder


## check chrmm project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT

if [[ ${#SLURM_ARRAY_TASK_ID} == 1 ]]
then 
	DIR="0${SLURM_ARRAY_TASK_ID}"
else
	DIR="${SLURM_ARRAY_TASK_ID}"
fi

echo "Loading config file"
source ./integrative/segway/config/config.txt
echo

#-----------------------------------------------------------------------#

module load Anaconda3
source activate segway

echo 'Changing to directory' ${LOADDIR}
cd ${LOADDIR}

FILES=$(find $LOADDIR -name "N+*.gnmdat")

if [ $# == 1 ] || [[ $2 =~ 'TRAIN' ]]
then
	mkdir -p ${TRAINDIR}

	cd ${SCRIPTDIR}
	sh ./integrative/segway/processing/learnModel.sh ${SLURM_ARRAY_TASK_ID} ${FILES}
fi


if [ $# == 1 ] || [[ $2 =~ 'ANNOTATE' ]]
then
	mkdir -p ${ANNODIR}
	mkdir -p ${AGGREDIR}

	cd ${SCRIPTDIR}
	sh ./integrative/segway/processing/annotateModel.sh ${SLURM_ARRAY_TASK_ID} ${FILES}
fi


end=`date +%s`

runtime=$((end-start))
echo 'Runtime:' $runtime

## move log files into a folder
cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}