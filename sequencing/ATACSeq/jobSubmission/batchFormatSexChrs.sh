#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/formatSexChr-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/formatSexChr-%A_%a.e
#SBATCH --job-name=formatSexChr-%A_%a.e

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source $1  

## reformat bam file
module load GATK
module load SAMtools

## load sample to process from text file

sampleName=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/Stage1Samples.txt | tail -1))

sh ./ATACSeq/preprocessing/11_subsetSexChrs.sh ${sampleName}

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv formatSexChr-${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}/