#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=general/logFiles/%u/download-%A_%a.e # error file
#SBATCH --output=general/logFiles/%u/download-%A_%a.o # output file
#SBATCH --job-name=download

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder


## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./sequencing/BSSeq/config/config.txt 
echo "Project directory is: " $DATADIR

## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#


mkdir -p ${RAWDATADIR}

link=$(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/dataLinks.txt | tail -n -1 | tr  -d '\r') # strip carriage return

echo ${link}


cd $RAWDATADIR
wget --no-check-certificate ${link}

## move log files into a folder
cd ${SCRIPTDIR}/general/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}