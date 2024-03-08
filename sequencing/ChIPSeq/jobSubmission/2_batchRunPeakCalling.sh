#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ChIPSeq/logFiles/%u/ChIPPeakCalling-%A_%a.o
#SBATCH --error=ChIPSeq/logFiles/%u/ChIPPeakCalling-%A_%a.e
#SBATCH --job-name=ChIPPeakCalling

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ChIPSeq/config/config.txt 
echo "Project directory is: " $DATADIR

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#

echo "Changing Folder to Data directory "
echo ${ALIGNEDDIR}

cd ${ALIGNEDDIR}

sampleName=($(head -n `expr ${SLURM_ARRAY_TASK_ID} + 1` ${METADIR}/samples.txt | tail -1))

if [[ $2 == 'control' ]]
then 
	column=$(awk -v RS=',' '/controlID/{print NR; exit}' $METADIR/sampleSheet.csv)
	control=$(awk -F"," -v col="$column" /$sampleName/'{print $col}' ${METADIR}/sampleSheet.csv)
	echo 'Control file specified:' $control
	shift
fi

if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then
	module purge
	module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
	module load BEDTools

	mkdir -p ${PEAKDIR}

	# derive histone mark from sampleSheet file
	## get column number for mark
	column=$(awk -v RS=',' '/target/{print NR; exit}' $METADIR/sampleSheet.csv)
	if [[ $column == '' ]]
	then
		{ echo 'target not found in sampleSheet, please check required columns'; exit 1; }
	fi

	echo 'Column number is ' $column
	## get mark for sample
	mark=$(cat ${METADIR}/sampleSheet.csv | awk -F"," -v col="$column" /$sampleName/'{print $col}')
	echo 'Mark is' $mark

	cd ${SCRIPTDIR}/
	sh ./ChIPSeq/preprocessing/samplePeaks.sh ${mark} ${sampleName} ${control}
	
	if [[ $? == 0 ]]
		then echo "Peaks called"
	fi
	
fi

if [ $# = 1 ] || [[ $2 =~ 'FRIP' ]]
then
	echo Starting calculate FRIP:
	date -u
	module purge
	module load BEDTools
	module load SAMtools

	mkdir -p ${PEAKDIR}/QCOutput

	cd ${SCRIPTDIR}/
	sh ./ChIPSeq/preprocessing/calcFrip.sh ${sampleName}

	if [[ $? == 0 ]]
		then date -u
		echo "FRIP calculated"
	fi
	
fi

echo "EXITCODE: " $?

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
