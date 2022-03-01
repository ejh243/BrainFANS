#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/ATACPeakCalling-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/ATACPeakCalling-%A_%a.e
#SBATCH --job-name=ATACPeakCalling-%A_%a.e

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source $1 

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
	./ATACSeq/preprocessing/5_shiftAlignedReads.sh ${sampleName}

	date -u
	echo "Reads shifted"
fi

if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then
	echo Starting peak calling at:
	date -u
	module purge
	module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
	module load BEDTools
	cd ${SCRIPTDIR}/
	./ATACSeq/preprocessing/6_samplePeaks.sh ${sampleName}
	
	date -u
	echo "Peaks called"
	
fi

if [ $# = 1 ] || [[ $2 =~ 'FRIP' ]]
then
	echo Starting calculate FRIP:
	date -u
	module purge
	module load BEDTools
	module load SAMtools

	sh ./ATACSeq/preprocessing/7_calcFrip.sh ${sampleName}
	date -u
	echo "FRIP calculated called"
	
fi


## move log files into a folder
mkdir -p LogFiles/ATAC/${SLURM_ARRAY_JOB_ID}
mv LogFiles/ATAC/ATACPeakCalling-${SLURM_ARRAY_JOB_ID}* LogFiles/ATAC/${SLURM_ARRAY_JOB_ID}