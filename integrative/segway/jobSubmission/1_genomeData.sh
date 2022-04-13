#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --mem=150G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/gnmdatGenerate-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/gnmdatGenerate-%A_%a.e
#SBATCH --job-name=gnmdatGenerate

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
source ./integrative/segway/config/config.txt

module load Anaconda3
source activate segway

sampleName=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/samples.txt | tail -1))
toProcess=($(find ${PEAKDIR} -maxdepth 1 -name ${sampleName}'*.filt'))

FILES=($(find ${PEAKDIR} -maxdepth 1 -name '*.filt'))

echo "File found is: " $( basename ${toProcess} )

sh ./integrative/segway/processing/bamtobed.sh ${sampleName} ${toProcess}