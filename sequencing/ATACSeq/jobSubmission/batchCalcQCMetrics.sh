#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/calcATACQC-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/calcATACQC-%A_%a.e
#SBATCH --job-name=calcATACQC-%A_%a

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1

module load R/3.6.3-foss-2020a

Rscript ${SCRIPTDIR}/ATACSeq/preprocessing/4_fragmentDistribution.r ${ALIGNEDDIR} ${SLURM_ARRAY_TASK_ID} 

## move log files into a folder
mkdir -p ATACSeq/logFiles/${SLURM_ARRAY_JOB_ID}
mv ATACSeq/logFiles/calcATACQC-${SLURM_ARRAY_JOB_ID}* ATACSeq/logFiles/${SLURM_ARRAY_JOB_ID}

