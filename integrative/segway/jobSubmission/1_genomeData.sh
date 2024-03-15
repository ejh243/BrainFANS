#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/gnmdatGenerate-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/gnmdatGenerate-%A_%a.e
#SBATCH --job-name=gnmdatGenerate

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder


## check chrmm project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT
PROJECT=$2
tissue=$3

#change

echo "Loading config file"
source ./integrative/segway/config/config.txt

#check tissue specified and assign default if not
if [[ $3 == '' ]]
then
	echo 'Tissue not specified, using default of prefrontal cortex|PFC'
	tissue="prefrontal cortex|PFC"
fi
echo

#if test -d ${PEAKDIR}
#then
#	echo 'Peak dir not found, using methyldir as ''...

#-----------------------------------------------------------------------#

module load Anaconda3
source activate segway

echo 'Changing to directory' $PEAKDIR
cd ${PEAKDIR}

## get the sample name and file to process
if test -f ${METADIR}/stage1Samples.txt 
then
	sampleName=($(head -n `expr ${SLURM_ARRAY_TASK_ID} + 1` ${METADIR}/stage1Samples.txt | tail -1))
	toProcess=($(find ${PEAKDIR} -name ${sampleName}'*.filt'))

	echo "File found is: " $( basename ${toProcess[0]} )
else
	echo 'stage1Samples.txt not found, exiting'
	exit 1
fi

## dont process samples outside of specified tissue
## get column number of tissue 
tissueCol=$(awk '$1 == "tissue"{print NR;exit} ' RS="," ${METADIR}/sampleSheet.csv)
echo 'Tissue column number is ' $tissueCol

sampleTissue=$(grep $sampleName ${METADIR}/sampleSheet.csv | awk -F',' -v col=$tissueCol '{print $col}')
if [[ $sampleTissue =~ $tissue ]]
then 
	echo 'Sample is in specified tissue:' $tissue
else
	echo "Sample is not in specified tissue, exiting." 

	cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
	mkdir -p ${SLURM_ARRAY_JOB_ID}
	mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}

	exit 1
fi

cd ${SCRIPTDIR}
sh ./integrative/segway/processing/genomedataLoad.sh ${sampleName} ${toProcess[0]}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}