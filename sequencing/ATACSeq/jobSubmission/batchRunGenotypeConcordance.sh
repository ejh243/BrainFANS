#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=LogFiles/ATAC/%u/verifyBAMID-%A_%a.o
#SBATCH --error=LogFiles/ATAC/%u/verifyBAMID-%A_%a.e
#SBATCH --job-name=verifyBAMID-%A_%a.e

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

## reformat bam file

module load picard/2.6.0-Java-1.8.0_131
module load GATK
module load SAMtools

# process a line from IDMap file
IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/matchedVCFIDs.txt | tail -1))

IDS=($(head -n 22 ${METADIR}/matchedVCFIDs.txt | tail -1))

sh ./ATACSeq/preprocessing/compareBamWithGenotypes.sh ${IDS[@]}


## move log files into a folder
mkdir -p logFiles/ATAC/$USER/${SLURM_ARRAY_JOB_ID}
mv logFiles/ATAC/$USER/verifyBAMID-${SLURM_ARRAY_JOB_ID}* logFiles/ATAC/$USER/${SLURM_ARRAY_JOB_ID}