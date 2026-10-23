#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPADownAnalysisS10-%A.log
#SBATCH --error=PIPPADownAnalysisS10-%A.err
#SBATCH --job-name=PIPPADownAnalysisS10

## ===================================================================================================================##
##                                      PIPPA STEP 10: Preparation for downstream analysis                            ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/10_batchDownstreamAnalyses.sh <project directory> <STEP> <GROUP> <THRES>        ||
##                                                                                                                    ||
## DESCRIPTION: This script prepares curated peak sets generated for each cell-type for downstream analyses.          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: PEAKS, COUNTS. Can be combined. Default is to run all                            ||
## $3 -> <GROUP> Specify cell group to run previous step(s).                                                          ||
## $4 -> <THRES> Specify peak calling threshold used for curated peak set.                                            ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/downAnalysisSamples_${GROUP}.txt that lists sample names for each cell type group.           ||
## - Config.txt, config.r files in <project directory>.                                                               ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PEAKCOUNTS ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA                                                       ||
## - A list of cell types of samples specified in config files: CELLTYPES                                             ||
## - A conda environment setup with several modules: samtools, MACS3, R                                               ||
## - For other tools, modules or references required, please refer to each subscript run in this script.              ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: curatedPeaksChecks.Rmd, countsInCleanPeaksCT.r, countsInCuratedPeaks.r                             ||
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
mv PIPPADownAnalysisS10-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

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

export GROUP=$3
export THRES=$4

cat << EOF
Available cell fractions are: ${CELLTYPES[@]}

EOF

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "SUMMARY" ]] && [[ ! $2 == "CTCLEAN" ]] && [[ ! $2 == 'COUNTS' ]];
then 
    { echo "Unknown step specified. Please use CTCLEAN, COUNTS or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## check cell type group has been specified for specific steps
if [[ ! " ${CELLTYPES[*]} " == *" $GROUP "* && ($2 == "COUNTS" )]]; 
then
  { echo "Cell group specified not found. Please check whether samples that belong to that cell fraction exist."; exit 1; }            
fi

## check peak calling threshold has been specified for specific steps
if [[ $2 == "COUNTS" && -z "${THRES:-}" ]];
then 
    { echo "Invalid threshold for $2 STEP. Should be a peak calling threshold used to generated a curated peak set from STEP 9.5" ; exit 1; }
fi

if [[ ! ${GROUP} == '' ]]
then
  echo "Chosen cell fraction is ${GROUP}"
 
  ## Retrieve samples from the same cell fraction that passed QC
  mapfile -t SAMPLES < "${META_DIR}/downAnalysisSamples_${GROUP}.txt"
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Number of samples found in ${GROUP} is:" """${#SAMPLES[@]}"""
  echo "Samples that belong to ${GROUP} are:" ${SAMPLES[@]}
else
  echo "No cell group specified."
fi

## ========= ##
##   STEPS   ##
## ========= ##

## option SUMMARY: summarise properties of curated peak sets generated at STEP 9.5 for all cell types
if [ $# = 1 ] || [[ $2 == 'SUMMARY' ]]
then

cat <<EOF

|| Running STEP 10.1 of PIPPA: SUMMARY ||

Properties of curated peak sets generated across cell types are checked.
Output directory is ${PEAK_DIR}/QCOutput

EOF

  Rscript -e "rmarkdown::render('${RSCRIPTS_DIR}/curatedPeaksChecks.Rmd', output_file='${PEAK_DIR}/QCOutput/curatedPeaksChecks.html')" "${CONFIGR}"

fi

## option CTCLEAN: counts in peaks for samples that pass all QC stages across cell types 
if [ $# = 1 ] || [[ $2 == 'CTCLEAN' ]]
then
  mkdir -p ${COUNTS_DIR_GROUPS}
  
cat <<EOF

|| Running STEP 10.2 of PIPPA: CTCLEAN ||

Counts in peak sets for all cell types will be performed.
Only samples that passed every QC step will be used.
Peaks generated for peaks curation will be used (STEP 9.3)
Output directory is ${COUNTS_DIR_GROUPS}

EOF

  Rscript "${RSCRIPTS_DIR}/countsInCleanPeaksCT.r" ${CONFIGR}
fi

## option COUNTS: counts in curated set of peaks for specified cell type
if [ $# = 1 ] || [[ $2 == 'COUNTS' ]]
then

  mkdir -p ${COUNTS_DIR_DOWN}
 
cat <<EOF

|| Running STEP 10.3 of ATAC-seq pipeline: COUNTS ||

Counts in curated peak sets on ${GROUP} will be performed.
Peak set with cutoff ${THRES} is chosen
Output directory is ${COUNTS_DIR_DOWN}

EOF

 Rscript "${RSCRIPTS_DIR}/countsInCuratedPeaks.r" ${CONFIGR} ${GROUP} ${THRES}
fi

echo Job finished on:
date -u

conda deactivate