#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ChIPSeq/logFiles/%u/calcChIPQC-%A.o
#SBATCH --error=ChIPSeq/logFiles/%u/calcChIPQC-%A.e
#SBATCH --job-name=calcChIPQC

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
JOBNAME="calcChIPQC"
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR


echo "Loading config file for project: " $1
export PROJECT=$1
source ./ChIPSeq/config/config.txt 

echo "Project directory is: " $DATADIR

#-----------------------------------------------------------------------#

mkdir -p ${PEAKDIR}/QCOutput

module load R/3.6.3-foss-2020a

Rscript ${SCRIPTDIR}/ChIPSeq/preprocessing/5_calcQCMetrics.r ${PROJECT}

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_JOB_ID}
mv ${JOBNAME}-${SLURM_JOB_ID}* ${SLURM_JOB_ID}