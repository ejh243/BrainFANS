#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSexChrS5-%A_%a.log
#SBATCH --error=ATACSexChrS5-%A_%a.err
#SBATCH --job-name=ATACSexChrS5

## ===================================================================================================================##
##                                 ATAC-seq pipeline STEP 5: Sex chromosomes                                          ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/5_batchFormatSexChrs.sh <project-directory> <STEP>        ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script prepares files for calling peaks in sex chromosomes, performs peak calling and checks the ||
## sex assigned to samples based on the sex chromosomes peaks                                                         ||
##                                                                                                                    || 
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: SPLIT, PEAKS, CHECK. Can be combined. Default is to run all                      ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, PEAKCOUNTS, ALIGNED_DIR, PEAK_DIR ||
##   SCRIPTS_DIR, PROJECT,RSCRIPTS_DIR                                                                                ||
## - Version/directory of the following modules should be specified in config file: SAMTVERS, RVERS, PIP_ENV          ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: subsetSexChrs.sh, sexChrPeaks.sh                                                                     ||
## - R scripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                                  ||
## - R Subscripts: collateSexChecks.r                                                                                 ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

if [[ $1 == '' ]] || [[ ! -d $1 ]]
then
  { echo "No project directory specified or could not be found." ; exit 1; }
else
  source "${1}/config.txt" 
fi

## load config file provided on command line related to the specified project
echo "Loading config file for project: " ${PROJECT}
echo "Project directory is: " $MAIN_DIR 
echo 'Script is running from directory: ' ${SCRIPTS_DIR}

LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACSexChrS5-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

## ================ ##
##    VARIABLES     ##
## ================ ##

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "SPLIT" ]] && [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "CHECK" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use SPLIT, PEAKS, CHECK or some combination of this as a single string (i.e. SPLIT,CHECK)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## If only $1 is specified, all steps will be run
## option SPLIT: subset reads from chromosomes X and Y to prepare samples for peak calling   
if [ $# = 1 ] || [[ $2 =~ 'SPLIT' ]]
then
  module purge
  module load $SAMTVERS
  
  echo " "
  echo "|| Running STEP 5.1 of ATAC-seq pipeline: SPLIT. Reads in sex chromosomes will be isolated.||"
  echo " "
  echo "Output directory will be: ${ALIGNED_DIR}/sexChr"
  echo " "
  
  ## load sample to process from samples that passed QC stage 1
  mapfile -t SAMPLEIDS < ${META_DIR}/samples.txt
  echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
  
  ## Sample(s) specified in by array number(s) from command line
  sampleID=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
  
  echo "Sample(s) specified by array number in command line is/are: " ${sampleID[@]}
  
  mkdir -p ${ALIGNED_DIR}/sexChr
  
  sh "${SUB_SCRIPTS_DIR}/subsetSexChrs.sh" ${sampleID}
  
fi

## option PEAKS: peaks are called only in sex chromosomes
if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then
  
  module purge
  module load $PVERS
  source ${PIP_ENV}/bin/activate
  module load $BEDTVERS
  
  echo " "
  echo "|| Running STEP 5.2 of ATAC-seq pipeline: PEAKS. Peaks will be called in sex chromosomes using MACS3 Single-end mode.||"
  echo " "
  echo "Output directory will be: ${PEAK_DIR}/ShiftedTagAlign/sexChr for peaks"
  echo "Output directory will be: ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/ for counts in peaks"
  echo " "
   
  mkdir -p ${PEAK_DIR}/ShiftedTagAlign/sexChr
  mkdir -p ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/
  
  
  sh "${SUB_SCRIPTS_DIR}/sexChrPeaks.sh"

fi

## option CHECK: results from previous steps are used to check assigned sex of samples
if [ $# = 1 ] || [[ $2 =~ 'CHECK' ]]
then
  
  module purge
  module load $RVERS
  
  echo " "
  echo "|| Running STEP 5.3 of ATAC-seq pipeline: CHECK. Sex of samples will be compared to be predicted sex.||"
  echo " "
  echo "Output directory will be: ${ALIGNED_DIR}/sexChr"
  echo " "
  
  Rscript ${RSCRIPTS_DIR}/collateSexChecks.r ${CONFIGR}
  
fi

echo Job finished on:
date -u
