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
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/jobSubmission/1_batchRunPreAnalysis.sh <config file> <option>      ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <config file> directory of config file for project analysed                                                  ||
## $2 -> <option> Specify step to run: FASTQC, TRIM, ALIGN or ENCODE. Can be combined. Default is to                  ||
## run all                                                                                                            ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the first step in the ATAC-seq pipeline: pre-analysis. This                      ||
## includes quality control of raw reads (FASTQC), trimming the reads (TRIM), aligning reads to                       ||
## reference genome (ALIGN) and calculate ENCODEQC metrics (ENCODE).                                                  ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${METADIR}/samples.txt that lists sample names.                                                          ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time when script runs
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source $1
echo "Loading config file for project: " ${PROJECT}

## check directories
echo "Project directory is: " $DATADIR 
echo "Script is running from directory: " ${SCRIPTDIR}

## Log files directory
LOG_DIR=ATACSeq/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACpreAnalysisS1-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
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
if test -f ${METADIR}/samples.txt;
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 
    echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
    
    ## Sample(s) specified in by array number(s) from command line
    sampleID=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
    toProcess=($(find ${RAWDATADIR} -maxdepth 1 -name ${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}'*'))
    
    ## sort the toProcess array so that R1 and R2 are consecutive 
    IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
    toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
    unset IFS 
    echo ${toProcess}
    echo "R1 file found is: " $( basename ${toProcess[0]} ) 
    echo "Sample to be used is: " ${sampleID} 
    
    ## If only $1 is specified, all steps will be run
    ## option FASTQC: run sequencing QC on raw reads files        
    if [ $# == 1 ] || [[ $2 =~ 'FASTQC' ]]
    then      
        echo "Running step 1.1 of ATAC-seq pipeline: FASTQC"
        module load FastQC 
        cd ${SCRIPTDIR}
        sh ./preScripts/fastqc.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}  
    fi
    
    ## option TRIM: run trimming on raw reads files and output trimmed reads.
    if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
    then
        module purge
        module load fastp
        
    	  echo "Running step 1.2 of ATAC-seq pipeline: TRIM"
        cd ${SCRIPTDIR}
        sh ./preScripts/fastp.sh ${sampleID} ${toProcess[0]} ${toProcess[1]} 
    fi

    ## option ALIGN: run alignment and post processing on sample
	  if [ $# == 1 ] || [[ $2 =~ 'ALIGN' ]]
  	then
   
  		echo "Running step 1.3 of ATAC-seq pipeline: ALIGN"
  		module purge ## had conflict issues if this wasn't run first
  		module load Bowtie2/2.3.4.2-foss-2018b
  		module load SAMtools
  		module load picard/2.6.0-Java-1.8.0_131
  		export PATH="$PATH:/lustre/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"
  		
  		cd ${SCRIPTDIR}
  		sh ./ATACSeq/preprocessing/alignment.sh ${sampleID}
  	fi
   
    ## option ENCODE: calculate ENCODEQC metrics on aligned samples
  	if [ $# == 1 ] || [[ $2 =~ 'ENCODE' ]]
  	then
      echo "Running step 1.4 of ATAC-seq pipeline: ENCODE"
  		module purge
  		module load SAMtools
  		module load BEDTools/2.27.1-foss-2018b ##necessary to specify earlier BEDTools version
          
      ## load conda env for samstats
      module load Anaconda3
      source $CONDAENV
      conda activate $ENVDIR

      cd ${SCRIPTDIR}
      sh ./ATACSeq/preprocessing/calcENCODEQCMetrics.sh ${sampleID}
    fi

else
    echo "File list not found"
fi


