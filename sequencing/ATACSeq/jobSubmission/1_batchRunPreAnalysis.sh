#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACpreAnalysisS1-%A_%a.log
#SBATCH --error=ATACpreAnalysisS1-%A_%a.err
#SBATCH --job-name=ATACpreAnalysisS1

## ===================================================================================================================##
##                               ATAC-seq pipeline STEP 1: Pre-analysis                                               ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/1_batchRunPreAnalysis.sh <project directory> <STEP>       ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the first step in the ATAC-seq pipeline: pre-analysis. This                      ||
## includes quality control of raw reads (FASTQC), trimming the reads (TRIM), aligning reads to                       ||
## reference genome (ALIGN) and calculate ENCODEQC metrics (ENCODE).                                                  ||
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
##   SCRIPTS_DIR, PROJECT                                                                                             ||
## - Version/directory of the following modules should be specified in config file: BOWTIEVERS, SAMTVERS, PICARDVERS, ||
##   BEDTVERS, ACVERS, CONDA_ENV, CONDA, UTILS.                                                                       ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: fastqc.sh, fastp.sh, alignment.sh and calcENCODEQCMetrics.sh                                         ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time when script runs
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source "${1}/config.txt"
echo "Loading config file for project: " ${PROJECT}

## check directories
echo "Project directory is: " $MAIN_DIR 
echo "Script is running from directory: " ${SCRIPTS_DIR}

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACpreAnalysisS1-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR


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
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] && [[ ! $2 =~ "ENCODE" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## check if file containing list of sample IDs exists and if so:
if test -f ${META_DIR}/samples.txt;
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < ${META_DIR}/samples.txt 
    echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
    
    ## Sample(s) specified in by array number(s) from command line
    sampleID=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
    toProcess=$(find ${RAWDATADIR} -maxdepth 1 -name ${sampleID}'*')
    ## sort the toProcess array so that R1 and R2 are consecutive 
    IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
    toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
    unset IFS 
    
    echo "R1 file found is: " $(basename ${toProcess[0]} )
    
    echo "Current sample: " ${sampleID} 
    
    ## option FASTQC: run sequencing QC on raw reads files        
    if [ $# == 1 ] || [[ $2 =~ 'FASTQC' ]]
    then
      module purge      
      module load FastQC 
      
      echo "|| Running STEP 1.1 of ATAC-seq pipeline: FASTQC ||"
      echo "Starting fastqc on ${sampleID} at: "
      date -u
      echo " "
      echo "Output written to " ${FASTQCDIR}
      
      sh "${SUB_SCRIPTS_DIR}/fastqc.sh" ${sampleID} ${toProcess[0]} ${toProcess[1]}  
          
    fi
    
    ## option TRIM: run trimming on raw reads files and output trimmed reads.
    if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
    then
      module purge
      module load fastp
      module load FastQC
      
  	  echo "|| Running STEP 1.2 of ATAC-seq pipeline: TRIM ||"
      echo "Starting trimming on ${sampleID} at: "
      date -u
      echo " "
      echo "Output written to " ${TRIM_DIR}
      
      sh "${SUB_SCRIPTS_DIR}/fastp.sh" ${sampleID} ${toProcess[0]} ${toProcess[1]}
        
    fi

    ## option ALIGN: run alignment and post processing on sample
	  if [ $# == 1 ] || [[ $2 =~ 'ALIGN' ]]
  	then
   
  		module purge
  		module load $BOWTIEVERS
  		module load $SAMTVERS
  		module load $PICARDVERS
  		export PATH="$PATH:${UTILS}"
     
      echo " "
      echo "|| Running STEP 1.3 of ATAC-seq pipeline: ALIGN ||"
      echo " "
      echo "Output written to ${ALIGNED_DIR} and ${ALIGNED_DIR}/QCOutput"
      
  		sh "${SUB_SCRIPTS_DIR}/alignment.sh" ${sampleID}
  	fi
   
    ## option ENCODE: calculate ENCODEQC metrics on aligned samples
  	if [ $# == 1 ] || [[ $2 =~ 'ENCODE' ]]
  	then
      
  		module purge
  		module load $BEDTVERS
      module load $SAMTVERS 
          
      ## load conda env for samstats
      module load ${ACVERS}
      source ${CONDA_ENV}
      conda activate ${CONDA}
      
      echo " "
      echo "|| Running STEP 1.4 of ATAC-seq pipeline: ENCODE ||"
      echo " "
      echo "Output written to ${ALIGNED_DIR}/ENCODEMetrics"
      
      sh "${SUB_SCRIPTS_DIR}/calcENCODEQCMetrics.sh" ${sampleID}
    fi

else
    echo "File list not found"
fi

echo Job finished on:
date -u 
