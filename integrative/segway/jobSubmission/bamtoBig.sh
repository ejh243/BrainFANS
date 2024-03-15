#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/bamtoBig-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/bamtoBig-%A_%a.e
#SBATCH --job-name=bamtoBig

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder


## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
DATA=$2
source ./sequencing/${DATA}/config/config.txt

mkdir -p ${ALIGNEDDIR}/coverage

module load deepTools

sampleName=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/samples.txt | tail -1))
toProcess=($(find ${ALIGNEDDIR} -maxdepth 1 -name ${sampleName}'*sorted.bam'))

echo "File found is: " $( basename ${toProcess} )

sh ./integrative/segway/processing/bamToBigwig.sh ${sampleName} ${toProcess}

## move log files into a folder
mkdir -p integrative/segway/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
mv integrative/segway/logFiles/${USER}/bamtoBig-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* integrative/segway/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
