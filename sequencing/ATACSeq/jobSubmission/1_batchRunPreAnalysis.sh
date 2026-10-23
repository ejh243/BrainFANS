#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPApreAnalysisS1-%A_%a.log
#SBATCH --error=PIPPApreAnalysisS1-%A_%a.err
#SBATCH --job-name=PIPPApreAnalysisS1

## ===================================================================================================================##
##                                           PIPPA STEP 1: Pre-analysis                                               ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/1_batchRunPreAnalysis.sh <project directory> <STEP>       ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the first step of PIPPA: pre-analysis. This  includes quality control of raw     ||
## reads (FASTQC), trimming the reads (TRIM), aligning reads to reference genome (ALIGN) and calculate ENCODEQC       ||
## metrics (ENCODE).                                                                                                  ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: FASTQC, TRIM, ALIGN or ENCODE. Can be combined. Default is to run all            ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, RAWDATADIR, ALIGNED_DIR, TRIM_DIR ||
##   SCRIPTS_DIR, PROJECT, CONDA_ENV, CONDA,  SAMPLES_LIST                                                            ||
## - A conda environment setup with several modules: samtools, fastqc, samstats, bedtools, bowtie2, picard and fastp  ||
## - For more software or references required, please refer to each subscript run in this script.                     ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: fastqc.sh, fastp.sh, alignment.sh, assign_multimappers.py and calcENCODEQCMetrics.sh                 ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time when script runs
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv PIPPApreAnalysisS1-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

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

## check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step has been specified, run all otherwise
if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "FASTQC" ]] && [[ ! $2 == "TRIM" ]] && [[ ! $2 == "ALIGN" ]] && [[ ! $2 == "ENCODE" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## check if file containing list of sample IDs exists and if so:
if test -f "${SAMPLES_LIST}";
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < "${SAMPLES_LIST}" 
    echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
    
    ## Sample(s) specified in by array number(s) from command line
    sampleName=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
    echo $sampleName
    toProcess=$(find ${RAWDATADIR} -maxdepth 1 -name ${sampleName}'*')
    
    ## sort the toProcess array so that R1 and R2 are consecutive 
    IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
    toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
    echo ${toProcess[@]}
    unset IFS 
    
    echo "R1 file found is: " $(basename ${toProcess[0]} )
    
    echo "Current sample: " ${sampleName} 
    
    ## option FASTQC: run sequencing QC on raw reads files        
    if [ $# == 1 ] || [[ $2 == 'FASTQC' ]]
    then
      
cat <<EOF

|| Running STEP 1.1 of PIPPA: FASTQC ||

Starting fastqc on ${sampleName}
Output written to ${FASTQCDIR}

EOF
      
      sh "${SUB_SCRIPTS_DIR}/fastqc.sh" ${sampleName} ${toProcess[0]} ${toProcess[1]}  
          
    fi
    
    ## option TRIM: run trimming on raw reads files and output trimmed reads.
    if [ $# == 1 ] || [[ $2 == 'TRIM' ]]
    then
      
cat <<EOF

|| Running STEP 1.2 of PIPPA: TRIM ||

Starting trimming on ${sampleName} 
Output written to ${TRIM_DIR}

EOF
      
      sh "${SUB_SCRIPTS_DIR}/fastp.sh" ${sampleName} ${toProcess[0]} ${toProcess[1]}
      
    fi

    ## option ALIGN: run alignment and post processing on sample
	  if [ $# == 1 ] || [[ $2 == 'ALIGN' ]]
  	then
   
    mkdir -p  ${ALIGNED_DIR}/QCOutput
	
cat <<EOF

|| Running STEP 1.3 of PIPPA: ALIGN ||

Starting trimming on ${sampleName}
Output written to ${ALIGNED_DIR} and ${ALIGNED_DIR}/QCOutput

EOF
      
      sh "${SUB_SCRIPTS_DIR}/alignment.sh" ${sampleName}
     
  	fi
   
    ## option ENCODE: calculate ENCODEQC metrics on aligned samples
  	if [ $# == 1 ] || [[ $2 == 'ENCODE' ]]
  	then
   
      mkdir -p ${ALIGNED_DIR}/ENCODEMetrics
      
cat <<EOF

|| Running STEP 1.4 of PIPPA: ENCODE ||

Output written to ${ALIGNED_DIR}/ENCODEMetrics

EOF
      
      sh "${SUB_SCRIPTS_DIR}/calcENCODEQCMetrics.sh" ${sampleName}

    fi

else
    echo "Sample list not found"
fi

conda deactivate 

echo Job finished on:
date -u 