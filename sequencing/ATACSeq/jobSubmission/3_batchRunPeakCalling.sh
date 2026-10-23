#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPAPeakCallingS3-%A_%a.log
#SBATCH --error=PIPPAPeakCallingS3-%A_%a.err
#SBATCH --job-name=PIPPAPeakCallingS3

## ===================================================================================================================##
##                                        PIPPA STEP 3: Sample Peak Calling                                           ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/3_batchRunPeakCalling.sh <project directory> <STEP>       ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the third step of PIPPA: calling peaks at sample level using the paired-end mode ||
## in MACS3 (PEAKS). Then results are collated to calculate Fraction Reads in Peaks (FRIP)                            ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input. Only for PEAKS. Set  ||
## to any number for FRIP.                                                                                            ||
## $1 -> <project directory> directory to config file for the corresponding project                                   ||
## $2 -> <STEP> Specify step to run: PEAKS, FRIP. Can be combined. Default is to run all. FRIP step should be run at  ||
## the end                                                                                                            ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PROJECT    ||
##   SCRIPTS_DIR, PEAK_DIR_SAMPLES, RSCRIPTS_DIR, CONDA_ENV, CONDA                                                    ||
## - For software or references required, please refer to each subscript run in this script.                          ||
## - A conda environment setup with several modules: samtools, MACS3, R, bedtools, samtools                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                            ||
## - Subscripts: shiftAlignedReads.sh, samplePeaks.sh, collateCalcFrip.sh                                             ||
##                                                                                                                    ||
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
mv PIPPAPeakCallingS3-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir:  $LOG_DIR
Running scripts from: $PWD

EOF

## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}

## ================ ##
##    VARIABLES     ##
## ================ ##

## check array specified and exit if not for BW and SCORE steps
if [[ $2 == "PEAKS" && ! -v SLURM_ARRAY_TASK_ID ]];
then 
    { echo "Invalid -array number specified for $2 STEP. Each number corresponds to a sample from the sample list provided." ; exit 1; }
fi

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "PEAKS" ]] && [[ ! $2 == "FRIP" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use PEAKS, FRIP or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## Collate all filtered, no duplicated, bam files
cd ${ALIGNED_DIR}
BAMFILES=($(ls *.filt.nodup.bam))
echo "Number of aligned bam files: "" ""${#BAMFILES[@]}"""

## Collate all samples
mapfile -t SAMPLEIDS < "${SAMPLES_LIST}"
echo "Total number of samples: "" ""${#SAMPLEIDS[@]}"""

## ============ ##
##    STEPS     ##
## ============ ##

## option PEAKS: peak calling is performed using MACS3 in paired-end at sample level in narrow peak mode 
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then

  mkdir -p ${PEAK_DIR_SAMPLES}
  
  ##Select sampleName using --array number
  sampleName=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
  echo "sampleName selected is ${sampleName}"

cat <<EOF

|| Running STEP 3.1 of PIPPA: PEAKS. Peaks will be called on sampleName using MACS3 on paired-end and narrowpeak mode. ||

sampleName specified by array number in command line is: ${sampleName}
Output directory is ${PEAK_DIR_SAMPLES}

EOF

  sh "${SUB_SCRIPTS_DIR}/samplePeaks.sh" ${sampleName}
fi

## option FRIP: collate peak calling results statistics in a single csv file
if [ $# = 1 ] || [[ $2 == 'FRIP' ]]
then
  
  mkdir -p ${PEAK_DIR}/QCOutput
  
cat <<EOF

|| Running STEP 3.2 of PIPPA: FRIP. Results from peak calling will be collated on a single file ||

Output directory is ${PEAK_DIR}/QCOutput
Samples found: ${SAMPLEIDS[@]}
Number samples found: ${#SAMPLEIDS[@]}

EOF

  sh "${SUB_SCRIPTS_DIR}/collateCalcFrip.sh" 

fi

conda deactivate
 
echo Job finished on:
date -u