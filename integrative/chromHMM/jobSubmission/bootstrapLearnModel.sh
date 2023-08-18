#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/bootLearn-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/bootLearn-%A_%a.e
#SBATCH --job-name=bootLearn

## print start date and time
echo Job started on:
date -u

start=`date +%s`


## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT


itNo=${SLURM_ARRAY_TASK_ID} 

echo 'Merged data directory is' ${MERGEDIR} 
echo 'Model directory is' ${MODELDIR}
echo

#-----------------------------------------------------------------------#

source ./integrative/chromHMM/config/config.txt

if [ $# == 1 ] || [[ $2 =~ 'MERGE' ]]
then
	module purge
	module load Java 

	mkdir -p ${MERGEDIR}

	echo "Running MergeBinary"

	java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar MergeBinary ${BINARISEDIR} ${MERGEDIR}

	echo "Finished binarising data"

fi

if [ $# == 1 ] || [[ $2 =~ 'MODEL' ]]
then
	module purge
	module load Java

	echo "Started learning model at:"
	date -u

	for stateNo in {5..15}
	do
		java -mx6000M -jar ${CHROMHMM}/ChromHMM.jar LearnModel -p 0 ${MERGEDIR} ${MODELDIR} ${stateNo} hg38
	done

	# note the job id for the likelihood stat later on
	echo ${SLURM_ARRAY_JOB_ID} >> ${MODELDIR}/jobid.txt

	if [[ $? == 0 ]]
	then 
		echo "Model(s) learned"
	fi
fi


echo 'EXIT CODE: ' $?

end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}
