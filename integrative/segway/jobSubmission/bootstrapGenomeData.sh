#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/bootLoad%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/bootLoad-%A_%a.e
#SBATCH --job-name=bootLoad

## print start date and time
echo Job started on:
date -u

start=`date +%s`


## check integrated project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT

export GROUP=N+

itNo=${SLURM_ARRAY_TASK_ID}
echo 'Loading config file for iteration number:' ${itNo}
source ./integrative/segway/config/config.txt
echo

#-----------------------------------------------------------------------#

module load Anaconda3
source activate segway


FILES=$(find ${LOADDIR} -name '*.bg')

for file in $FILES 
do
	## get prefix for 
	sampleName=$(basename $(sed 's/[.].*//' <<< $file))

	echo 'PREFIX is:' ${sampleName}

	filename=$( basename $file )

	cd ${SCRIPTDIR}
	sh ./integrative/segway/processing/genomedataLoad.sh ${sampleName} ${filename}
done



## move log files into a folder
cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}