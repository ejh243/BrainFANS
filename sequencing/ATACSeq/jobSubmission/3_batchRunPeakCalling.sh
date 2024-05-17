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
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/jobSubmission/3_batchRunPeakCalling.sh <project directory> <option>||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                       ||
## $1 -> <project directory> directory to config file for the corresponding project                                   ||
## $2 -> <option> Specify step to run: SHIFT, PEAKS, FRIP. Can be combined. Default is to run all                     ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the core analysis of the ATAC-seq pipeline, which is calling peaks at sample     ||
## level using the paired-end mode in MACS3.                                                                          ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source "${1}/config.txt"
echo "Loading config file for project: " $1
echo "Project directory is: " $DATADIR

## Log files directory
LOG_DIR=ATACSeq/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACPeakCallingS3-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

## check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "SHIFT" ]] && [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "FRIP" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use SHIFT, PEAKS, FRIP or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## Collate all filtered, no duplicated, bam files
cd ${ALIGNEDDIR}
BAMFILES=($(ls *.filt.nodup.bam))
echo "Number of aligned bam files: "" ""${#BAMFILES[@]}"""
  
## ============ ##
##    STEPS     ##
## ============ ##

## If only $1 is specified, all steps will be run
## option SHIFT: shift aligned reads for single-end peak calling   
if [ $# = 1 ] || [[ $2 =~ 'SHIFT' ]]
then
  module purge
	module load BEDTools/2.29.2-GCC-9.3.0
  module load SAMtools/1.11-GCC-9.3.0
	module load R/3.6.3-foss-2020a
 
	echo "Step 3.1 SHIFT started. Aligned reads will be shifted"
  mapfile -t SAMPLES < ${METADIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample(s) specified by array number in command line is/are: " ${sample}
  
  if [[ -s ${ALIGNEDDIR}/${sample}.filt.nodup.bam ]]
  then 
    cd ${SCRIPTDIR}
	  sh ./ATACSeq/preprocessing/shiftAlignedReads.sh ${sample}
  else
    echo "Aligned bam file for ${sample} not found. Aligned reads can't be shifted."
  fi
fi

## option PEAKS: peak calling is performed using MACS3 in paired-end 
if [ $# = 1 ] || [[ $2 =~ 'PEAKS' ]]
then

	module purge
	module load ${PVERS}
  source ${PIPENV}/bin/activate
	module load BEDTools
 
  echo "Step 3.2 PEAKS started. Aligned reads will be used for peak calling using MACS3 PE and TA"
  mapfile -t SAMPLES < ${METADIR}/samples.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample(s) specified by array number in command line is/are: " ${sample}
  
  mkdir -p ${PEAKDIR}/MACS/BAMPE
  
  if [[ -s ${ALIGNEDDIR}/${sample}.filt.nodup.bam ]]
  then 
    cd ${SCRIPTDIR}
    sh ./ATACSeq/preprocessing/samplePeaks.sh ${sample}
	else
    echo "Aligned bam file for ${sample} not found. Peaks can't be called in aligned reads."
  fi
fi

## option FRIP: collate peak calling results statistics in a single csv file
if [ $# = 1 ] || [[ $2 =~ 'FRIP' ]]
then

	module purge
	module load BEDTools/2.29.2-GCC-9.3.0
  module load SAMtools/1.11-GCC-9.3.0
  
  echo "Step 3.3 FRIP started. Peak calling statistics will be collated."
  mapfile -t SAMPLES < ${METADIR}/samples.txt
  
	mkdir -p ${PEAKDIR}/QCOutput

	cd ${SCRIPTDIR}/

  ## Outputs a single csv file with all peak calling results for all samples
  sh ./ATACSeq/preprocessing/collateCalcFrip.sh ${SAMPLES[@]}
	
fi

