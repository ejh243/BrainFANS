#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPAcalcQCS2-%A_%a.log
#SBATCH --error=PIPPAcalcQCS2-%A_%a.err
#SBATCH --job-name=PIPPAcalcQCS2

## ===========================================================================================================================##
##                                         PIPPA STEP 2: Post alignment processing                                            ##
## ===========================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/2_batchCalcQCMetrics.sh <project directory> <STEP> <GROUP_SAMPS>  ||
##                                                                                                                            ||
## DESCRIPTION: This script performs the second step of PIPPA: ATAC-seq-specific metrics calculation. Includes the            ||
## classification of aligned reads into fragment sizes and related statistics (FRAG).Then pre-processing of aligned reads are || 
## done (BW) to calculate TSS enrichment score (SCORE) and collate results (COLLATE)                                          ||
##                                                                                                                            ||
## INPUTS:                                                                                                                    || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input. Needed for steps BW,         ||
##  and SCORE. Set to any number for FRAG or COLLATE                                                                          || 
## - $1 -> <project directory> : directory to config file for project                                                         ||
## - $2 -> <STEP> Specify step to run: FRAG, BW, SCORE or COLLATE. Needs to be specified.                                     ||
## - $3 -> <GROUP_SAMPS> Only to be specified for FRAG step. It splits the samples into groups of 10 to run in                ||
##  parallel, so batch number should be the number of samples divided by 10. E.G. If there are 50 samples, batch              ||
##  number should be 0-5, producing 5 batches of 10 samples each. If the number of samples is less than 10, set batch         ||
##  number to 0.                                                                                                              ||
##                                                                                                                            ||
## REQUIRES:                                                                                                                  ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                                 ||
## - Config.txt file in <project directory>.                                                                                  ||
## - The following variables specified in config file: MAIN_DIR, SCRIPTS_DIR, LOG_DIR, ALIGNED_DIR, RSCRIPTS_DIR,             ||
##   CONDA, CONDA_ENV                                                                                                         ||
## - For other variables, software or references required, please refer to each subscript run in this script.                 ||
## - A conda environment setup with several modules: samtools, R, samstats,bedtools, bowtie2 and fastp                        ||
## - fragmentDistribution.r file in ${RSCRIPTS_DIR} = ./Rscripts directory                                                    ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                                    ||
## - Subscripts: prepareBigWig.sh, calcTssEnrich.sh, collateTSSScores.sh                                                      ||
##                                                                                                                            ||
## ===========================================================================================================================##

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
mv PIPPAcalcQCS2-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

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

export $BATCH_NUM=$3

if [[ $2 == '' ]] 
then
  { echo "No step specified. All steps in the script will be run. "; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 == "FRAG" ]] && [[ ! $2 == "BW" ]] && [[ ! $2 == "SCORE" ]] && [[ ! $2 == "COLLATE" ]] &&  [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FRAG, BW, SCORE, COLLATE" ; exit 1; }            
fi

## check array specified and exit if not for BW and SCORE steps
if [[ ( $2 == "FRAG" || $2 == "BW" || $2 == "SCORE" ) && ! -v SLURM_ARRAY_TASK_ID ]];
then 
    { echo "Invalid -array number specified for $2 STEP. Each number corresponds to a sample from the sample list provided." ; exit 1; }
fi

## check group of samples specified and exit if not for FRAG step
if [[ $2 == "FRAG" && -z "${BATCH_NUM:-}" ]]
then
  { echo "No group of samples specified for STEP 2.1 FRAG. "; exit 1; }
fi

## Collate all samples
mapfile -t SAMPLESIDS < "${SAMPLES_LIST}"

echo "Total number of samples: "" ""${#SAMPLESIDS[@]}"""

##  ==================  ##
##     CALCQCMetrics    ##
##  ==================  ##

## option FRAG: Calculate fragment distribution of batch of samples 
if [ $# == 1 ] || [[ $2 == 'FRAG' ]]
then

  mkdir -p ${ALIGNED_DIR}/QCOutput

cat <<EOF

|| Running STEP 2.1 of PIPPA: Post alignment processing (QC metrics and fragment distribution). ||

Calculating QC metrics and fragment distribution for samples in group of samples ${BATCH_NUM} 
Output directory is ${ALIGNED_DIR}/QCOutput/

EOF
  
  Rscript "${RSCRIPTS_DIR}/fragmentDistribution.r" ${CONFIGR} ${BATCH_NUM} 
  
  if [[ ! -f "${ALIGNED_DIR}/QCOutput/FSD_${BATCH_NUM}.pdf" ]] && [[ ! -f "${ALIGNED_DIR}/QCOutput/FragmentDistribution_${BATCH_NUM}.rdata" ]]
  then 
    { echo "Fragment distribution and other QC metrics of samples in batch ${BATCH_NUM} could not be calculated. Please make sure STEP 1.3 ALIGN was successfully run."  ; exit 1;}
  else
    echo "Fragment distribution and other QC metrics of samples in batch ${BATCH_NUM} were successfully calculated."
  fi
fi

## option BW: Convert  bam file to bigWig format
if [ $# == 1 ] || [[ $2 == 'BW' ]]
then

  export TSS_DIR="${ALIGNED_DIR}/bigWig/tssEnrich/"
  mkdir -p "${TSS_DIR}"
  
  ##Select sample using --array number
  sampleName=${SAMPLESIDS[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample selected is ${sampleName}"

cat <<EOF

|| Running STEP 2.2 of PIPPA: Prepare file for TSS calculation. ||

Transforming sample ${sampleName} to bigWig format to calculate TSS enrichment score. 
Output directory is ${TSS_DIR}

EOF

  sh "${SUB_SCRIPTS_DIR}/prepareBigWig.sh" ${sampleName}
  
fi

## option SCORE: Calculate TSS enrichment score
if [ $# == 1 ] || [[ $2 == 'SCORE' ]]
then

  mkdir -p "${TSS_DIR}"
  
  ##Select sampleName using --array number
  sampleName=${SAMPLESIDS[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample selected is ${sampleName}"

cat <<EOF

|| Running STEP 2.3 of PIPPA: Calculate TSS enrichment. ||

Calculate TSS enrichment score for sample ${sampleName}.
Output directory is ${TSS_DIR}

EOF

  sh "${SUB_SCRIPTS_DIR}/calcTssEnrich.sh" ${sampleName}

fi

## option COLLATE: Collate TSS enrichment scores for all samples
if [ $# == 1 ] || [[ $2 == 'COLLATE' ]]
then

  mkdir -p ${ALIGNED_DIR}/QCOutput  
  
cat <<EOF

|| Running STEP 2.4 of PIPPA: Collate TSS enrichment scores. ||

TSS enrichment scores for all samples are collated in a single file.
Output directory is ${ALIGNED_DIR}/QCOutput 

EOF

  sh "${SUB_SCRIPTS_DIR}/collateTSSscores.sh" 

fi

conda deactivate
 
echo Job finished on:
date -u