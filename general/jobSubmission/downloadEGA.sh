#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/downloadEGA-%A_%a.o
#SBATCH --error=general/logFiles/downloadEGA-%A_%a.e
#SBATCH --job-name=downloadEGA

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
SUBPROJECT=$2

source ./sequencing/BSSeq/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}
echo "metadata directory is: " $METADIR

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#=============================================================================#
module purge
module load Python/3.9.6-GCCcore-11.2.0
source ../.ebi/bin/activate


CONFIG=${SCRIPTDIR}/general/config/egaCredentials.json

downloadIDs=($(grep "fastq" ${METADIR}/${SUBPROJECT}/delimited_maps/Sample_File.map | awk -F'\t' '{print $4}'))

echo "Number of download ids is:" ${#downloadIDs[@]}

echo "to process is:" ${downloadIDs[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${RAWDATADIR}
cd ${RAWDATADIR}

echo "Downloading files"

pyega3 -cf ${CONFIG} fetch ${downloadIDs[$SLURM_ARRAY_TASK_ID]}