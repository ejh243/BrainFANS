#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/calcATACQC-%A_%a.o
#SBATCH --error=LogFiles/calcATACQC-%A_%a.e
#SBATCH --job-name=calcATACQC-%A_%a
#SBATCH --array=0-29%10 

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

module load R/3.6.3-foss-2020a

Rscript ${SCRIPTDIR}/ATACSeq/batchCalcATACSeqQCMetrics.r ${SCRIPTDIR}/ATACSeq/config.r ${SLURM_ARRAY_TASK_ID} 