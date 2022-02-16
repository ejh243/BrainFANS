#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/ATACQCSummary-%A.o
#SBATCH --error=ATACSeq/logFiles/ATACQCSummary-%A.e
#SBATCH --job-name=ATACQCSummary


## print start date and time
echo Job started on:
date -u


## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

module load MultiQC
## use multiqc to collate QC output statistics
'''
mkdir -p ${FASTQCDIR}/multiqc
cd ${FASTQCDIR}/
multiqc . -f -o ${FASTQCDIR}/multiqc

mkdir -p ${ALIGNEDDIR}/multiqc
cd ${ALIGNEDDIR}/
multiqc . -f -o ${ALIGNEDDIR}/multiqc
'''
## run other bespoke utilty scripts to collate other QC metrics
cd ${SCRIPTDIR}/

#./ATACSeq/preprocessing/8_progressReport.sh 
./ATACSeq/preprocessing/9_countMTReads.sh 
#./ATACSeq/preprocessing/10_collateFlagStatOutput.sh 

## move log files into a folder
mkdir -p ATACSeq/logFiles/${USER}
mv ATACSeq/logFiles/ATACQCSummary-${SLURM_ARRAY_JOB_ID}* ATACSeq/logFiles/${USER}
