#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=96:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=integrative/chromTools/logFiles/%u/complete-%A_%a.e # error file
#SBATCH --output=integrative/chromTools/logFiles/%u/complete-%A_%a.o # output file
#SBATCH --job-name=complete

#-----------------------------------------------------------------------#

## needs to be run with array number 0-number of targets 
## at the moment can only be run from $USER=jms260 !!! (hard coded chromtools path)

## print start date and time
echo Job started on:
date -u

PROJECT=$1
source ./integrative/chromTools/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#

## column of targets
mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 

## convert to BED file

echo "sampleids are" ${SAMPLEIDS[@]}
sampleid=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}

module purge
module load Python/3.7.2-GCCcore-6.4.0
source ~/.chromTools/bin/activate

files=${sampleid}.bed

echo "Files are: " $files

mkdir -p $TEMPDIR

CHROMTOOLDIR=${CHROMTOOLDIR}/replicate/$sampleid

mkdir -p $CHROMTOOLDIR

cd $ALIGNEDDIR

chromTools complete -f $files -o ${CHROMTOOLDIR} --force-overwrite -i 5_000_000


## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromTools/logFiles/$USER
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}