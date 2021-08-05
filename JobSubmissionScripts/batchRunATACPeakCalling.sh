#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/ATAQPeakCalling-%A_%a.o
#SBATCH --error=LogFiles/ATAQPeakCalling-%A_%a.e
#SBATCH --job-name=ATAQPeakCalling-%A_%a.e
#SBATCH --array=0-295%40 ## runs 19 jobs with 40 at any one time


## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

echo "Changing Folder to Data directory "
echo ${ALIGNEDDIR}

cd ${ALIGNEDDIR}
BAMFILES=($(ls *_depDup_q30.bam))

echo "Number of bam files found for alignment:"" ""${#BAMFILES[@]}"""	

sample=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%_depDup_q30.bam})

## shift reads for peak calling
echo "Shifting reads"

module load BEDTools
module load SAMtools
module load R/3.6.3-foss-2020a

cd ${SCRIPTDIR}
./ATACSeq/shiftAlignedReadsPE.sh ${sampleName}

date -u
echo "Reads shifted"


echo Starting peak calling at:
date -u
module purge
module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
module load BEDTools
cd ${SCRIPTDIR}/
./ATACSeq/peakCallingPE.sh ${sampleName}
