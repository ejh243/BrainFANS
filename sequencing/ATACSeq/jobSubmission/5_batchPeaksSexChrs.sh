#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPASexChrS5-%A_%a.log
#SBATCH --error=PIPPASexChrS5-%A_%a.err
#SBATCH --job-name=PIPPASexChrS5

## ===================================================================================================================##
##                                         PIPPA STEP 5: Sex chromosomes                                              ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/5_batchFormatSexChrs.sh <project-directory> <STEP>        ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the fifth step of PIPPA and first step of Stage 2 of QC. It prepares files for   ||
## calling peaks in sex chromosomes (SPLIT) and performs peak calling (PEAKS) on sex chromosome reads.                ||
##                                                                                                                    || 
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input for SPLIT. Set to any ||
## number of PEAKS                                                                                                    || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: SPLIT, PEAKS                                                                     ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR_SEXCHR,     ||
##   SCRIPTS_DIR, PROJECT,RSCRIPTS_DIR, CONDA, CONDA_ENV, PEAKCOUNTS_DIR_SEXCHR                                       ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, MACS3                                                  ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: subsetSexChrs.sh, sexChrPeaks.sh                                                                     ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv PIPPASexChrS5-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

## load config file provided on command line related to the specified project
cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR 
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir: $LOG_DIR
Running scripts from: $PWD

EOF

## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}

## ================ ##
##    VARIABLES     ##
## ================ ##

##check array specified and exit if not
if [[ $2 == "SPLIT" && ! -v SLURM_ARRAY_TASK_ID ]];
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "SPLIT" ]] && [[ ! $2 == "PEAKS" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use SPLIT, PEAKS or some combination of this as a single string (i.e. SPLIT, PEAKS)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## If only $1 is specified, all steps will be run
## option SPLIT: subset reads from chromosomes X and Y to prepare samples for peak calling   
if [ $# = 1 ] || [[ $2 == 'SPLIT' ]]
then
  
  ## load sample to process from samples that passed QC stage 1
  mapfile -t SAMPLEIDS < "${SAMPLES_LIST}"
  ## Sample(s) specified in by array number(s) from command line
  sampleName=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}

cat <<EOF

|| Running STEP 5.1 of PIPPA: SPLIT. Reads in sex chromosomes will be isolated.||

Output directory will be: ${ALIGNED_DIR}/sexChr
Number of sample IDs found: ${#SAMPLEIDS[@]}
Sample specified by array number in command line is:  ${sampleName[@]}

EOF
  
  mkdir -p ${ALIGNED_DIR}/sexChr
  
  sh "${SUB_SCRIPTS_DIR}/subsetSexChrs.sh" ${sampleName}

fi

## option PEAKS: peaks are called only in sex chromosomes
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
	
cat <<EOF

|| Running STEP 5.2 of PIPPA: PEAKS. Peaks will be called in sex chromosomes using MACS3 Paired-end mode.||

Output directory will be: ${PEAK_DIR_SEXCHR} for peaks in sex chromosomes
Output directory will be: ${COUNTS_DIR_SEXCHR} for counts in peaks

EOF

  mapfile -t SAMPLEIDS < "${SAMPLES_LIST}"
  echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""

  mkdir -p ${PEAK_DIR_SEXCHR}
  mkdir -p ${COUNTS_DIR_SEXCHR}

  sh "${SUB_SCRIPTS_DIR}/sexChrPeaks.sh"
  
fi

conda deactivate 

echo Job finished on:
date -u