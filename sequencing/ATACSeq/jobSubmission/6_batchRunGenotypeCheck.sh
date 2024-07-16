#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACGenoConcorS6-%A_%a.o
#SBATCH --error=ATACGenoConcorS6-%A_%a.e
#SBATCH --job-name=ATACGenoConcorS6

## ============================================================================================================================##
##                              ATAC-seq pipeline STEP 6: Genotype concordance                                                 ##
## ============================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/6_batchRunGenotypeConcordance.sh <project directory> <STEP>        ||
## - execute from pipeline's main directory                                                                                    ||
##                                                                                                                             ||
## DESCRIPTION: This script will detect possible DNA contamination in order to ensure high quality sequence reads. If any      ||
## contamination is detected, possible swaps will be suggested.                                                                ||
##                                                                                                                             ||
## INPUTS:                                                                                                                     || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input                                ||
## $1 -> <project directory> Directory of the project, location of the config file                                             ||
## $2 -> <STEP> Specify step to run: COMPARE, GENCHECK. Can be combined. Default is to run all                                 ||
##                                                                                                                             ||
## OUTPUTS:                                                                                                                    ||
## - Potential contaminated samples will be specified in: ${META_DIR}/potentialSwitches.txt                                    ||
##                                                                                                                             ||
## REQUIRES:                                                                                                                   ||
## - File with samples IDs and their matching VCF IDs: ${META_DIR}/matchedVCFIDs.txt                                           ||
## - Config.txt file in <project directory>.                                                                                   ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, SCRIPTS_DIR, PROJECT,PEAK_DIR ||
## - Version/directory of the following modules should be specified in config file: PICARDVERS, RVERS, SAMTVERS,GATKVERS       ||
## - Softwares: Picard, GATK, SAMtools, R, Pandoc,                                                                             ||
## - For modules or references required, please refer to each subscript run in this script.                                    ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                                     ||
## - Subscripts: compareBamWithGenotypes.sh, searchBestGenoMatch.sh                                                            ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                                        ||
## - R scripts: collateSampleChecks.Rmd                                                                                        ||
##                                                                                                                             ||
## ============================================================================================================================##

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
mv ATACGenoConcorS6-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

## ================ ##
##    VARIABLES     ##
## ================ ##

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "COMPARE" ]] && [[ ! $2 =~ "GENCHECK" ]] && [[ ! $2 =~ "SWITCH" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use COMPARE, GENCHECK, SWITCH or some combination of this as a single string (i.e. SPLIT,CHECK)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## option COMPARE: genotype of samples are verified by comparing them with existing genotype information   
if [ $# = 1 ] || [[ $2 =~ 'COMPARE' ]]
then
  module purge
  module load $PICARDVERS
  module load $GATKVERS
  module load $SAMTVERS
  
  echo " "
  echo "|| Running STEP 6.1 of ATAC-seq pipeline: COMPARE. Samples will be compared to their matching genotype data.||"
  echo " "
  echo "Output directory will be: ${ALIGNED_DIR}/genotypeConcordance"
  echo "Output directory will be: ${ALIGNED_DIR}/baseRecalibrate"
  echo " "
  
  # process a line from IDMap file
  IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${META_DIR}/matchedVCFIDs.txt | tail -1))
  echo "Samples to verify with their corresponding VCF IDs are: ${IDS[@]}"
  
  sh "${SUB_SCRIPTS_DIR}/compareBamWithGenotypes.sh" ${IDS[@]}
  
fi

## option GENCHECK: Sex check and Genotype check results are collated in Rmarkdown   
if [ $# = 1 ] || [[ $2 =~ 'GENCHECK' ]]
then
  
	module purge
	module load ${RVERS}
	module load Pandoc
 
  echo " "
  echo "|| Running STEP 6.2 of ATAC-seq pipeline: GENCHECK. Results from genotype and sex check will be collated in a Rmarkdown report.||"
  echo " "
  echo "Output directory will be: ${PEAK_DIR}/QCOutput"
  echo " "
 
	Rscript -e "rmarkdown::render(paste0(commandArgs(trailingOnly=TRUE)[1], '/collateSampleChecks.Rmd'), output_file=paste0(commandArgs(trailingOnly=TRUE)[2], '/QCOutput/stage2SummaryStats.html'))" "${RSCRIPTS_DIR}" "$PEAK_DIR" "${CONFIGR}"
  
fi

## option SWITCH: Samples with non-matching genotype are suggested to be switched
if [ $# = 1 ] || [[ $2 =~ 'SWITCH' ]]
then
  module purge
  module load $PICARDVERS
  module load $GATKVERS
  module load $SAMTVERS
  
  echo " "
  echo "|| Running STEP 6.3 of ATAC-seq pipeline: SWITCH. Samples with potential genotype contamination will be selected for potential switches.||"
  echo " "
  echo "Output directory will be: ${ALIGNED_DIR}/genotypeConcordance"
  echo "Output directory will be: ${ALIGNED_DIR}/baseRecalibrate"
  echo " "
  
  cd ${ALIGNED_DIR}/genotypeConcordance/
  
  awk '{if($7 != "FREEMIX") print FILENAME,$0}' *.selfSM > collatedGenoCheckStats.txt
  awk '{if($7 != "FREEMIX" && ($12 > 0.9)) print FILENAME,$1}' *.selfSM > ${META_DIR}/potentialSwitches.txt

  IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${META_DIR}/potentialSwitches.txt | tail -1))
  echo "Samples that are likely to be contaminated are: ${ID[@]}" 
  
  sh "${SUB_SCRIPTS_DIR}/searchBestGenoMatch.sh" ${IDS[@]}
 
fi

echo Job finished on:
date -u
