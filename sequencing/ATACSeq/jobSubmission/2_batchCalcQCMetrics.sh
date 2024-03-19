#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACcalcQCS2-%A_%a.log
#SBATCH --error=ATACcalcQCS2-%A_%a.err
#SBATCH --job-name=ATACcalcQCS2

## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 2: Post alignment processing                                    ##
## ===================================================================================================================##
## PROOFREADER : Marina Flores-Payan                                                                                  ||
## CONTACT: m.flores-payan@exeter.ac.uk                                                                               ||
## LAST UPDATE: February 2024                                                                                         ||
## ===================================================================================================================##
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/jobSubmission/2_batchCalcQCMetrics.sh <config file>                ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> It splits the samples into groups of 10 to run in parallel, so batch number should be the number of     ||
##  samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples|| 
##  each. If the number of samples is less than 10, set batch number to 0                                             ||
## - $1 -> <config file> : directory to config file for project                                                       ||
##                                                                                                                    ||
## DESCRIPTION: This script uses the ATACseqQC R package to generate the fragment distribution, calculate periodicity ||
## and get ratios of nucleosomefree, mono, bi, tri reads.                                                             ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source $1
echo "Loading config file for project: " $1
echo "Project directory is: " $DATADIR

## Log files directory
LOG_DIR=ATACSeq/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv "ATACcalcQCS2-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}*" $LOG_DIR

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

##  ==================  ##
##     CALCQCMetrics    ##
##  ==================  ##

mkdir -p ${ALIGNEDDIR}/QCOutput

module load R/4.2.1-foss-2022a

echo "Running step 2 of ATAC-seq pipeline: Post alignment processing (QC metrics and fragment distribution)."
echo "Calculating QC metrics and fragment distribution for samples in batch " ${SLURM_ARRAY_TASK_ID} 
echo "Output directory is " "${ALIGNEDDIR}/QCOutput/"
Rscript ${SCRIPTDIR}/ATACSeq/preprocessing/fragmentDistribution.r $PROJECT ${SLURM_ARRAY_TASK_ID} 

echo "QC metrics and fragment distribution calculated"
