#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
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
##   CONDA, CONDA_ENV
## - Directory of the following software should be specified in config file: verifyBAMID (VERIFYBAMID)                         ||
## - Softwares: Picard, GATK, R, Pandoc, samtools in conda environment                                                         ||
## - For modules or references required, please refer to each subscript run in this script.                                    ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                                     ||
## - Subscripts: compareBamWithGenotypes.sh, searchBestGenoMatch.sh                                                            ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                                        ||
## - R scripts: collateSampleChecks.Rmd                                                                                        ||
## - --array number should start from 1 for the first sample.                                                                  ||
## - For STEP 6.2 GENCHECK, array number should be a single number, as need to be run once.                                    ||
##                                                                                                                             ||
## ============================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv ATACGenoConcorS6-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

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
if [[ ! $2 =~ "COMPARE" ]] && [[ ! $2 =~ "GENCHECK" ]] && [[ ! $2 =~ "SWITCH" ]] && [[ ! $2 =~ "COLLATE" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use COMPARE, COLLATE, GENCHECK, SWITCH or some combination of this as a single string (i.e. SPLIT,CHECK)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## option COMPARE: genotype of samples are verified by comparing them with existing genotype information   
if [ $# = 1 ] || [[ $2 =~ 'COMPARE' ]]
then
  
cat <<EOF

|| Running STEP 6.1 of ATAC-seq pipeline: COMPARE. ||

Samples will be compared to their matching genotype data.

Output directory will be: ${ALIGNED_DIR}/genotypeConcordance
Output directory will be: ${ALIGNED_DIR}/baseRecalibrate

EOF
  
  # process a line from IDMap file
  IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${META_DIR}/matchedVCFIDs.txt | tail -1))
  echo "Samples to verify with their corresponding VCF IDs are: ${IDS[@]}"
  
  sh "${SUB_SCRIPTS_DIR}/compareBamWithGenotypes.sh" ${IDS[@]}
  
fi

## option COLLATE: collate other QC metrics
if [ $# = 1 ] || [[ $2 =~ 'COLLATE' ]]
then

cat <<EOF

|| Running step 6.2 of ATAC-seq pipeline: COLLATE ||

It will be checked that all samples had gone through all processes and had produced the right outputs.
Output directory is ${META_DIR}

EOF
 
	sh ${SUB_SCRIPTS_DIR}/progressReportS2.sh 
	
fi


## option GENCHECK: Sex check and Genotype check results are collated in Rmarkdown   
if [ $# = 1 ] || [[ $2 =~ 'GENCHECK' ]]
then

cat <<EOF

|| Running STEP 6.3 of ATAC-seq pipeline: GENCHECK. ||

Results from genotype and sex check will be collated in a Rmarkdown report.
Output directory will be: ${PEAK_DIR}/QCOutput

EOF
 
	Rscript -e "rmarkdown::render('${RSCRIPTS_DIR}/collateSampleChecks.Rmd', output_file='$PEAK_DIR/QCOutput/stage2SummaryStats.html')" "${CONFIGR}"
  
fi

## option SWITCH: Samples with non-matching genotype are suggested to be switched
if [ $# = 1 ] || [[ $2 =~ 'SWITCH' ]]
then
  
cat <<EOF

|| Running STEP 6.4 of ATAC-seq pipeline: SWITCH. ||

Samples with potential genotype contamination will be selected for potential switches.
Output directory will be: ${ALIGNED_DIR}/genotypeConcordance
Output directory will be: ${ALIGNED_DIR}/baseRecalibrate

EOF
  
  cd ${ALIGNED_DIR}/genotypeConcordance/
  
  awk '{if($7 != "FREEMIX") print FILENAME,$0}' *.selfSM > collatedGenoCheckStats.txt
  awk '{if($7 != "FREEMIX" && ($12 > 0.9)) print FILENAME,$1}' *.selfSM > ${META_DIR}/potentialSwitches.txt

  mapfile -t SAMPLEIDS < ${META_DIR}/potentialSwitches.txt
  echo "Number of samples likely to have been swapped: ${#SAMPLEIDS[@]}"
  
  IDS=($(head -n ${SLURM_ARRAY_TASK_ID} ${META_DIR}/potentialSwitches.txt | tail -1))
  echo "Samples that is likely to have been swapped with another sample: ${IDS[@]}" 
  
  sh "${SUB_SCRIPTS_DIR}/searchBestGenoMatch.sh" ${IDS[@]}

fi

conda deactivate

echo Job finished on:
date -u
