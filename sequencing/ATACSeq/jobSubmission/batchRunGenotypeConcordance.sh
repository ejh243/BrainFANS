#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/verifyBAMID-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/verifyBAMID-%A_%a.e
#SBATCH --job-name=verifyBAMID-%A_%a.e

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

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR

## reformat bam file

module load picard/2.6.0-Java-1.8.0_131
module load GATK
module load SAMtools

# process a line from IDMap file
IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/matchedVCFIDs.txt | tail -1))

IDS=($(head -n 22 ${METADIR}/matchedVCFIDs.txt | tail -1))

sh ./ATACSeq/preprocessing/compareBamWithGenotypes.sh ${IDS[@]}

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv verifyBAMID-${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}/