#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/genoSearch-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/genoSearch-%A_%a.e
#SBATCH --job-name=genoSearch

#-----------------------------------------------------------------------#

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

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#

## take samples that look contaminated and search for best

cd ${ALIGNEDDIR}/genotypeConcordance/
awk '{if($7 != "FREEMIX") print FILENAME,$0}' *.selfSM > collatedGenoCheckStats.txt
awk '{if($7 != "FREEMIX" && ($5 > 1000000 && $12 > 0.9)) print FILENAME,$1}' *.selfSM > ${METADIR}/potentialSwitches.txt


IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/potentialSwitches.txt | tail -1))

cd ${SCRIPTDIR}

sh ./ATACSeq/preprocessing/searchBestGenoMatch.sh ${IDS[@]}

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}