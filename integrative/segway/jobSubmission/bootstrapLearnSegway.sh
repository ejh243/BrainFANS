#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/bootLearn-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/bootLearn-%A_%a.e
#SBATCH --job-name=bootLearn

## print start date and time
echo Job started on:
date -u

start=`date +%s`


## check chrmm project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT

##check states have been specified specified and exit if not
if [[ $# == 1 ]]
then 
    { echo "Max and min number of states have not been specified. Please give CL arguments 2 and 3 as min and max states, exiting." ; exit 1; }
fi

echo 'Min state is' $2
echo 'Max state is' $3

export GROUP=N+

itNo=${SLURM_ARRAY_TASK_ID}
echo 'Loading config file for iteration number:' ${itNo}
source ./integrative/segway/config/config.txt

states=$(seq $2 $3)
echo 'States are:' $states
echo
#-----------------------------------------------------------------------#

module load Anaconda3
source activate segway

FILES=$(find $LOADDIR -name $GROUP"*.gnmdat")

for state in ${states}
do
	echo
	echo "State is" $state
	if [[ ${#state} == 1 ]]
	then 
		DIR="0${state}"
	else
		DIR="${state}"
	fi

	source ./integrative/segway/config/config.txt

	mkdir -p ${TRAINDIR}

	cd ${SCRIPTDIR}
	#sh ./integrative/segway/processing/learnModel.sh ${state} ${FILES}


	mkdir -p ${ANNODIR}
	mkdir -p ${AGGREDIR}

	cd ${SCRIPTDIR}
	sh ./integrative/segway/processing/annotateModel.sh ${state} ${FILES}

	end=`date +%s`

	runtime=$((end-start))
	echo 'Runtime for' $state 'states is:' $runtime
done


end=`date +%s`

runtime=$((end-start))
echo 'Total runtime:' $runtime

## move log files into a folder
cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}