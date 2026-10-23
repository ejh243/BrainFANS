#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPAGenoConcorS6-%A_%a.log
#SBATCH --error=PIPPAGenoConcorS6-%A_%a.err
#SBATCH --job-name=PIPPAGenoConcorS6

## ============================================================================================================================##
##                                          PIPPA STEP 6: Genotype concordance                                                 ##
## ============================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/6_batchRunGenotypeConcordance.sh <project directory> <STEP>        ||
##                                                                                                                             ||
## DESCRIPTION: This script checks for possible DNA contamination or sample switch in order to ensure high quality sequence    ||
## reads (COMPARE). If any contamination is detected, possible swaps will be suggested.                                        ||
##                                                                                                                             ||
## INPUTS:                                                                                                                     || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input.                               ||
## $1 -> <project directory> Directory of the project, location of the config file                                             ||
## $2 -> <STEP> Specify step to run: COMPARE, SWITCH. Can be combined. Default is to run all                                   ||
##                                                                                                                             ||
## OUTPUTS:                                                                                                                    ||
## - Potential contaminated samples will be specified in: ${META_DIR}/potentialSwitches.txt                                    ||
##                                                                                                                             ||
## REQUIRES:                                                                                                                   ||
## - File with samples IDs and their matching VCF IDs: ${META_DIR}/matchedVCFIDs.txt                                           ||
## - Config.txt file in <project directory>.                                                                                   ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, SCRIPTS_DIR, PROJECT,PEAK_DIR ||
##   CONDA, CONDA_ENV                                                                                                          ||
## - Softwares: Picard, GATK, R, Pandoc, samtools, verifyBAMID in conda environment                                            ||
## - For modules or references required, please refer to each subscript run in this script.                                    ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                                     ||
## - Subscripts: compareBamWithGenotypes.sh, searchBestGenoMatch.sh                                                            ||
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
mv PIPPAGenoConcorS6-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

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
if [[ ! -v SLURM_ARRAY_TASK_ID ]];
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 == "COMPARE" ]] && [[ ! $2 == "SWITCH" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use COMPARE, SWITCH or some combination of this as a single string (i.e. SPLIT,CHECK)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## option COMPARE: genotype of samples are verified by comparing them with existing genotype information   
if [ $# = 1 ] || [[ $2 == 'COMPARE' ]]
then
  
cat <<EOF

|| Running STEP 6.1 of PIPPA: COMPARE. ||

Samples will be compared to their matching genotype data.

Output directory will be: ${ALIGNED_DIR}/genotypeConcordance
Output directory will be: ${ALIGNED_DIR}/baseRecalibrate

EOF

  ## If directories does not exist, create
  mkdir -p "${ALIGNED_DIR}/baseRecalibrate"
  mkdir -p "${ALIGNED_DIR}/genotypeConcordance"

  if [ ! -s "${META_DIR}/matchedVCFIDs.txt" ]; then
  awk -F',' 'BEGIN {OFS="\t"}
  NR==1 {
    for (i=1; i<=NF; i++) {
      col = $i
      gsub(/\r/,"",col)

      if (col=="sampleID") s=i
      if (col=="vcfID") v=i
    }
    if (!s || !v) {
      print "DEBUG HEADER:"
      for (i=1; i<=NF; i++) print "[" $i "]" > "/dev/stderr"
      print "ERROR: columns not found" > "/dev/stderr"
      exit 1
    }
    next
  }
  { print $s, $v }
  ' "${SAMPLE_SHEET}" > "${META_DIR}/matchedVCFIDs.txt"
  fi
  
  read -r sampleID vcfID < <(
    sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${META_DIR}/matchedVCFIDs.txt"
  )


  echo "Samples to verify with their corresponding VCF IDs are: $sampleID $vcfID"

  sh "${SUB_SCRIPTS_DIR}/compareBamWithGenotypes.sh" "$sampleID" "$vcfID"
  
fi

## option SWITCH: Samples with non-matching genotype are suggested to be switched
if [ $# = 1 ] || [[ $2 == 'SWITCH' ]]
then
  
cat <<EOF

|| Running STEP 6.2 of PIPPA: SWITCH. ||

Samples with potential genotype contamination will be selected for potential switches.
Output directory will be: ${ALIGNED_DIR}/genotypeConcordance
Output directory will be: ${ALIGNED_DIR}/baseRecalibrate

EOF
  
  cd ${ALIGNED_DIR}/genotypeConcordance/
  
  awk '{if($7 != "FREEMIX") print FILENAME,$0}' *.selfSM > collatedGenoCheckStats.txt
  awk '{if($7 != "FREEMIX" && ($12 > 0.9)) print FILENAME,$1}' *.selfSM > "${META_DIR}/potentialSwitches.txt"

  mapfile -t SAMPLEIDS < "${META_DIR}/potentialSwitches.txt"
  echo "Number of samples likely to have been swapped: ${#SAMPLEIDS[@]}"
  
  IDS=($(head -n ${SLURM_ARRAY_TASK_ID} "${META_DIR}/potentialSwitches.txt" | tail -1))
  echo "Samples that is likely to have been swapped with another sample: ${IDS[@]}" 
  
  sh "${SUB_SCRIPTS_DIR}/searchBestGenoMatch.sh" ${IDS[@]}

fi

conda deactivate

echo Job finished on:
date -u