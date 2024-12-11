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
## DESCRIPTION: This script performs the core analysis of the ATAC-seq pipeline, which is calling peaks at sample     ||
## level using the paired-end mode in MACS3.                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <project directory> directory to config file for the corresponding project                                   ||
## $2 -> <STEP> Specify step to run: SHIFT, PEAKS, FRIP. Can be combined. Default is to run all                       ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PROJECT    ||
##   SCRIPTS_DIR, PEAK_DIR, RSCRIPTS_DIR, CONDA_ENV, CONDA                                                            ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, MACS3, R, bedtools, samtools                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                            ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - Subscripts: shiftAlignedReads.sh, samplePeaks.sh, collateCalcFrip.sh                                             ||
## - For STEP 3.3 FRIP, a single job array number should be used, e.g. -array=0                                       ||
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
mv ATACPeakCallingS3-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir:  $LOG_DIR

EOF


## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}

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
 
  mapfile -t SAMPLES < ${META_DIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  
cat <<EOF

|| Running STEP 3.1 of ATAC-seq pipeline: SHIFT. Reads will be shifted for later peak calling on sex chromosomes ||

Output directory is ${ALIGNED_DIR}
Sample specified by array number in command line is: ${sample}

EOF
  
  sh "${SUB_SCRIPTS_DIR}/shiftAlignedReads.sh" ${sample}

fi

## option PEAKS: peak calling is performed using MACS3 in paired-end 
if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then
	
  mapfile -t SAMPLES < ${META_DIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  mkdir -p ${PEAK_DIR}/BAMPE
  
cat <<EOF

|| Running STEP 3.2 of ATAC-seq pipeline: PEAKS. Peaks will be called on sample using MACS3 on Paired-end mode. ||

Output directory is ${PEAK_DIR}/BAMPE
Sample specified by array number in command line is: ${sample}

EOF
  
  sh "${SUB_SCRIPTS_DIR}/samplePeaks.sh" ${sample}
  
fi

## option FRIP: collate peak calling results statistics in a single csv file
if [ $# = 1 ] || [[ $2 =~ 'FRIP' ]]
then
  
  mapfile -t SAMPLES < ${META_DIR}/samples.txt
  mkdir -p ${PEAK_DIR}/QCOutput
  
cat <<EOF

|| Running STEP 3.3 of ATAC-seq pipeline: FRIP. Results from peak calling will be collated on a single file ||

Output directory is ${PEAK_DIR}/QCOutput
Samples found: ${SAMPLES[@]}
Number samples found: ${#SAMPLES[@]}

EOF

  sh "${SUB_SCRIPTS_DIR}/collateCalcFrip.sh" ${SAMPLES[@]}
  

fi

conda deactivate
 
echo Job finished on:
date -u
