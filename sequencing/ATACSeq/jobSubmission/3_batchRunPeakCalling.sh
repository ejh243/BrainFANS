#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/ATACPeakCalling-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/ATACPeakCalling-%A_%a.e
#SBATCH --job-name=ATACPeakCalling

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
JOBNAME="ATACPeakCalling"
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
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
BAMFILES=($(ls *.filt.nodup.bam))

echo "Number of bam files found for alignment:"" ""${#BAMFILES[@]}"""	

sample=${BAMFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%.filt.nodup.bam})


if [ $# = 1 ] || [[ $2 =~ 'SHIFT' ]]
then
	module load BEDTools
	module load SAMtools
	module load R/3.6.3-foss-2020a
	## shift reads for peak calling
	echo "Shifting reads"

	cd ${SCRIPTDIR}
	./ATACSeq/preprocessing/4_shiftAlignedReads.sh ${sampleName}

	if [[ $? == 0 ]]
		then date -u
		echo "Reads shifted"
	fi
fi

if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then
	echo Starting peak calling at:
	date -u
	module purge
	module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
	module load BEDTools
	cd ${SCRIPTDIR}/
	./ATACSeq/preprocessing/5_samplePeaks.sh ${sampleName}
	
	if [[ $? == 0 ]]
		then date -u
		echo "Peaks called"
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
	sh ./ATACSeq/preprocessing/6_calcFrip.sh ${sampleName}

	if [[ $? == 0 ]]
		then date -u
		echo "FRIP calculated called"
	fi
	
fi

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv ${JOBNAME}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
