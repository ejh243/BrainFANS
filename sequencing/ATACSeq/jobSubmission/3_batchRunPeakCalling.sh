#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACPeakCallingS3-%A_%a.log
#SBATCH --error=ATACPeakCallingS3-%A_%a.err
#SBATCH --job-name=ATACPeakCallingS3

## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 3: Sample Peak Calling                                          ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/3_batchRunPeakCalling.sh <project directory> <STEP>       ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <project directory> directory to config file for the corresponding project                                   ||
## $2 -> <STEP> Specify step to run: SHIFT, PEAKS, FRIP. Can be combined. Default is to run all                       ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the core analysis of the ATAC-seq pipeline, which is calling peaks at sample     ||
## level using the paired-end mode in MACS3.                                                                          ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PROJECT    ||
##   SCRIPTS_DIR, PEAK_DIR, RSCRIPTS_DIR                                                                              ||
## - Version/directory of the following modules should be specified in config file: RVERS, SAMTVERS, PIP_ENV, PVERS   ||
##   BEDTVERS, PIP_ENV                                                                                                ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                            ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - Subscripts: shiftAlignedReads.sh, samplePeaks.sh, collateCalcFrip.sh                                             ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source "${1}/config.txt"
echo "Loading config file for project: " $PROJECT
echo "Project directory is: " $MAIN_DIR

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACPeakCallingS3-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

## ================ ##
##    VARIABLES     ##
## ================ ##

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Invalid -array number specified. Each number corresponds to a sample from the sample list provided." ; exit 1; }
fi

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "SHIFT" ]] && [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "FRIP" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use SHIFT, PEAKS, FRIP or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## Collate all filtered, no duplicated, bam files
cd ${ALIGNED_DIR}
BAMFILES=($(ls *.filt.nodup.bam))
echo "Number of aligned bam files: "" ""${#BAMFILES[@]}"""
  
## ============ ##
##    STEPS     ##
## ============ ##

## option SHIFT: shift aligned reads for single-end peak calling   
if [ $# = 1 ] || [[ $2 =~ 'SHIFT' ]]
then
  module purge
	module load $BEDTVERS
  module load $SAMTVERS
	module load $RVERS
 
	echo " "
  echo "|| Running STEP 3.1 of ATAC-seq pipeline: SHIFT. Reads will be shifted for later peak calling on sex chromosomes ||"
  echo " "
  echo "Output directory is ${ALIGNED_DIR}"
  
  mapfile -t SAMPLES < ${META_DIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample(s) specified by array number in command line is/are: " ${sample}
  
  sh "${SUB_SCRIPTS_DIR}/shiftAlignedReads.sh" ${sample}
  
fi

## option PEAKS: peak calling is performed using MACS3 in paired-end 
if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then

	module purge
	module load ${PVERS}
  source ${PIP_ENV}/bin/activate
	module load $BEDTVERS
 
  echo " "
  echo "|| Running STEP 3.2 of ATAC-seq pipeline: PEAKS. Peaks will be called on sample using MACS3 on Paired-end mode. ||"
  echo " "
  echo "Output directory is ${PEAK_DIR}"
  
  mapfile -t SAMPLES < ${META_DIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample(s) specified by array number in command line is/are: " ${sample}
  
  mkdir -p ${PEAK_DIR}/BAMPE
  
  sh "${SUB_SCRIPTS_DIR}/samplePeaks.sh" ${sample}
	
fi

## option FRIP: collate peak calling results statistics in a single csv file
if [ $# = 1 ] || [[ $2 =~ 'FRIP' ]]
then

	module purge
	module load $BEDTVERS
  module load $SAMTVERS
  
  echo " "
  echo "|| Running STEP 3.3 of ATAC-seq pipeline: FRIP. Results from peak calling will be collated on a single file ||"
  echo " "
  echo "Output directory is ${PEAK_DIR}/QCOutput"
  
  mapfile -t SAMPLES < ${META_DIR}/samplesBatch2.txt
  echo "There are ${#SAMPLES[@]} samples found."
  
	mkdir -p ${PEAK_DIR}/QCOutput

  #sh "${SUB_SCRIPTS_DIR}/collateCalcFrip.sh" ${SAMPLES[@]}
  sh "/lustre/home/mf662/scripts_ISCA/subscripts/collateCalcFrip.sh" ${SAMPLES[@]}
	
fi

echo Job finished on:
date -u
