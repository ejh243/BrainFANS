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
## EXECUTION: sbatch --array= ./jobSubmission/2_batchCalcQCMetrics.sh <project directory>                             ||
## - execute from pipeline's main scripts directory                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> It splits the samples into groups of 10 to run in parallel, so batch number should be the number of     ||
##  samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples|| 
##  each. If the number of samples is less than 10, set batch number to 0                                             ||
##  The samples in the batch specified will be taken from the aligned samples in ALIGNED_DIR                          ||
## - $1 -> <project directory> : directory to config file for project                                                 ||
##                                                                                                                    ||
## DESCRIPTION: This script uses the ATACseqQC R package to generate the fragment distribution and calculate some     ||
## summary statistics (e.g. multimodality of reads, periodicity)                                                      ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: MAIN_DIR, SCRIPTS_DIR, LOG_DIR, ALIGNED_DIR, RSCRIPTS_DIR, CONDA, CONDA_ENV            ||
## - R version should be > 4.2 in conda environment.                                                                  ||
## - fragmentDistribution.r file in ${RSCRIPTS_DIR} = ./Rscripts directory                                            ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv ATACcalcQCS2-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR


cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir:  $LOG_DIR

EOF

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi


## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}

##  ==================  ##
##     CALCQCMetrics    ##
##  ==================  ##

mkdir -p ${ALIGNED_DIR}/QCOutput

cat <<EOF

|| Running STEP 2 of ATAC-seq pipeline: Post alignment processing (QC metrics and fragment distribution). ||

Calculating QC metrics and fragment distribution for samples in batch ${SLURM_ARRAY_TASK_ID} 
Output directory is ${ALIGNED_DIR}/QCOutput/

EOF

Rscript ${RSCRIPTS_DIR}/fragmentDistribution.r ${CONFIGR} ${SLURM_ARRAY_TASK_ID} 

if [[ ! -f ${ALIGNED_DIR}/QCOutput/FSD_batch_${SLURM_ARRAY_TASK_ID}.pdf ]] && [[ ! -f ${ALIGNED_DIR}/QCOutput/FragmentDistribution_Batch_${SLURM_ARRAY_TASK_ID}.rdata ]]
then 
  { echo "Fragment distribution and other QC metrics of samples in batch ${SLURM_ARRAY_TASK_ID} could not be calculated. Please make sure STEP 1.3 ALIGN was successfully run."  ; exit 1;}
else
  echo "Fragment distribution and other QC metrics of samples in batch ${SLURM_ARRAY_TASK_ID} were successfully calculated."
fi
