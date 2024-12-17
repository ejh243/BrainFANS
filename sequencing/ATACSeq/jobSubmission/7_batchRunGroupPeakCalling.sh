#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACGroupPeakCallingS7-%A.log
#SBATCH --error=ATACGroupPeakCallingS7-%A.err
#SBATCH --job-name=ATACGroupPeakCallingS7

## ===================================================================================================================##
##                           ATAC-seq pipeline STEP 7: Group Peak calling                                             ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/7_batchRunGroupPeakCalling.sh <project directory> <STEP> <GROUP>                ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script performs peak calling by grouping samples of the same cell type. This includes selecting  ||
## samples that passed both stages 1 and 2 of quality control and analysing their fragment distribution to double     ||
## they are good quality ATAC-seq samples (FRAGSIZE). Next, peak calling is performed (PEAKS). Counts in this peaks   ||
## are obtained for advanced analysis (COUNTS)                                                                        ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: FRAGSIZE, PEAKS, COUNTS. Can be combined. Default is to run all                  ||
## $3 -> <GROUP> Specify cell group to run FRAGSIZE or PEAKS steps. Not needed for COUNTS step.                       ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PEAKCOUNTS ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA                                                       ||
## - A list of cell types of samples specified in config file: CELLTYPES                                              ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, MACS3                                                  ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: groupPeaks.sh                                                                                        ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: fragDistributionForPeaks.r, countsInPeaks.r, samplesForGroupAnalysis.r                             ||
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
mv ATACGroupPeakCallingS7-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

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

export GROUP=$3

cat << EOF
Available cell fractions are: ${CELLTYPES[@]}

EOF

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "FRAGSIZE" ]] && [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "COUNTS" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FRAGSIZE, PEAKS, COUNTS or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

if [[ ! " ${CELLTYPES[*]} " == *" $GROUP "* && ($2 =~ "FRAGSIZE" || $2 =~ "PEAKS")]]; 
then
  { echo "Cell group specified not found. Please check whether samples that belong to that cell fraction exist."; exit 1; }            
fi


if [[ ! ${GROUP} == '' && ($2 =~ "FRAGSIZE" || $2 =~ "PEAKS")]]
then
  echo "Chosen cell fraction is ${GROUP}"
  ## Samples for group peak calling are gathered for input cell fraction
  Rscript ${RSCRIPTS_DIR}/samplesForGroupAnalysis.r ${CONFIGR} ${GROUP}

  ## Retrieve samples from the same cell fraction
  mapfile -t SAMPLES < ${META_DIR}/samplesForGroupAnalysisOrdered_${GROUP}.txt

  echo "Number of samples found in ${GROUP} is:" """${#SAMPLES[@]}"""
  echo "Samples that belong to ${GROUP} are:" ${SAMPLES[@]}
else
  echo "No cell group specified."
fi


## ========= ##
##   STEPS   ##
## ========= ##

## option FRAGSIZE: plot fragment size distribution of samples chosen for group peak calling to check they are good ATAC-seq samples
if [ $# = 1 ] || [[ $2 == 'FRAGSIZE' ]]
then
  
cat <<EOF

|| Running STEP 7.1 of ATAC-seq pipeline: FRAGSIZE ||

Fragment size distribution of samples will be output in a pdf file
Output directory is ${PEAKCOUNTS}

EOF

  Rscript "${RSCRIPTS_DIR}/fragDistributionForPeaks.r" ${CONFIGR} ${GROUP}
  
fi

## option PEAKS: peak calling is performed using MACS3 in paired-end at group level  
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
  
  mkdir -p ${PEAK_DIR}/BAMPE/group
  
cat <<EOF

|| Running STEP 7.2 of ATAC-seq pipeline: PEAKS ||

Peak calling on $GROUP samples will be performed.
Output directory is ${PEAK_DIR}/BAMPE/group

EOF

  sh "${SUB_SCRIPTS_DIR}/groupPeaks.sh" ${GROUP}

fi

## option COUNTS: get read counts in peaks called above
if [ $# = 1 ] || [[ $2 == 'COUNTS' ]]
then
  
  mkdir -p ${PEAKCOUNTS}/Counts/
  
cat <<EOF

|| Running STEP 7.3 of ATAC-seq pipeline: COUNTS ||

Counts in peaks will be output in a single file.
Output directory is ${PEAKCOUNTS}/Counts/

EOF

  Rscript "${RSCRIPTS_DIR}/countsInPeaks.r" ${CONFIGR} 
fi

conda deactivate

echo Job finished on:
date -u
