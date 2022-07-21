#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/sexCheck-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/sexCheck-%A_%a.e
#SBATCH --job-name=sexCheck

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR


## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}

#-----------------------------------------------------------------------#

if [ ! -d ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr ]
then
	## call peaks for sex chromosomes & do read counts in these peaks

	module load MACS2
	module load BEDTools
	sh ./ATACSeq/preprocessing/sexChrPeaks.sh
fi

module purge
module load R/3.6.3-foss-2020a
Rscript ATACSeq/preprocessing/collateSexChecks.r ${DATADIR}/

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}