#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --mem=150G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/test-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/test-%A_%a.e
#SBATCH --job-name=test

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
source /lustre/home/jms260/BrainFANS/sequencing/ATACSeq/config/config.txt

module load Anaconda3
source activate segway

sampleName=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/samples.txt | tail -1))
toProcess=($(find ${RAWDATADIR} -maxdepth 1 -name ${sampleName}'*'))

## sort the toProcess array so that R1 and R2 are consecutive 
IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
unset IFS 

echo "Raw r1 file found is: " $( basename ${toProcess[0]} )

sh ./integrative/segway/processing/bamtobed.sh ${sampleName} ${toProcess[0]} ${toProcess[1]}