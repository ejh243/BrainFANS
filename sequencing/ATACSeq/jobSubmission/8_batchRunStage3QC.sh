#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPAGroupPeakCallingS8-%A.log
#SBATCH --error=PIPPAGroupPeakCallingS8-%A.err
#SBATCH --job-name=PIPPAGroupPeakCallingS8

## ===================================================================================================================##
##                                       PIPPA STEP 8: Group Peak calling                                             ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/8_batchRunStage3QC.sh <project directory> <STEP> <GROUP>                        ||
##                                                                                                                    ||
## DESCRIPTION: This script performs peak calling by grouping samples of the same cell type. This includes selecting  ||
## samples that passed both stages 1 and 2 of quality control. Next, peak calling is performed (PEAKS). Counts in     ||
## this peaks are obtained for advanced analysis (COUNTS). Samples are clustered based on their cell-type assigned to ||
## detect potential failed FANS samples (CTCHECK).                                                                    ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: MERGE, PEAKS, COUNTS and CTCHECK. Can be combined.                               ||
##  Default is to run all                                                                                             ||
## $3 -> <GROUP> Specify cell group to run MERGE or PEAKS steps. Not needed for other steps.                          ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR_GROUPS,     ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA, PEAKCOUNTS_DIR_GROUPS                                ||
## - A list of cell types of samples specified in config file: CELLTYPES                                              ||
## - A conda environment setup with several tools: samtools, MACS3                                                    ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: groupPeaks.sh                                                                                        ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: countsInPeaksGroup.r, samplesForGroupAnalysis.r ,collateCellTypeCheck.Rmd                          ||
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
mv PIPPAGroupPeakCallingS8-${SLURM_ARRAY_JOB_ID}* $LOG_DIR


cat <<EOF

Loading config file for project: ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory: ${SCRIPTS_DIR}
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

cat << EOF
Available cell fractions are: ${CELLTYPES[@]}

EOF

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "MERGE" ]] && [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "COUNTS" ]] && [[ ! $2 == 'CTCHECK' ]] && [[ ! $2 == ' ' ]];
then 
    { echo "Unknown step specified. Please use MERGE, PEAKS, COUNTS or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

if [[ ! " ${CELLTYPES[*]} " == *" $GROUP "* && ($2 =~ "MERGE" || $2 =~ "PEAKS" )]]; 
then
  { echo "Cell group specified not found. Please check whether samples that belong to that cell fraction exist."; exit 1; }            
fi


if [[ ! ${GROUP} == '' && ( $2 =~ "MERGE" || $2 =~ "PEAKS")]];
then
  echo "Chosen cell fraction is ${GROUP}"
  ## Samples for group peak calling are gathered for input cell fraction
  Rscript "${RSCRIPTS_DIR}/samplesForGroupAnalysis.r" ${CONFIGR} ${GROUP}

  ## Retrieve samples from the same cell fraction
  mapfile -t SAMPLES < "${META_DIR}/samplesForGroupAnalysis_${GROUP}.txt"
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Number of samples found in ${GROUP} is:" """${#SAMPLES[@]}"""
  echo "Samples that belong to ${GROUP} are:" ${SAMPLES[@]}
else
  echo "No cell group specified."
fi


## ========= ##
##   STEPS   ##
## ========= ##

## option MERGE: sample from the same cell-type are merged into a single bam file
if [ $# == 1 ] || [[ $2 == 'MERGE' ]]
then

  OUPUTDIR="${ALIGNED_DIR}/merged/passS12/"
  mkdir -p "${OUPUTDIR}"
  
cat <<EOF

|| Running STEP 8.1 of PIPPA: MERGE ||

A merged bam file will be created with all samples for the chosen cell-type.
Output directory is ${OUPUTDIR}
Aligned files taken from ${ALIGNED_DIR}

EOF
  bamFiles=()
  for sampleName in ${SAMPLES[@]}
  do
    sampleBAM="${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"
    bamFiles+=( ${sampleBAM})
  done
  
  echo "Number of aligned bam files found in ${ALIGNED_DIR} is: ${#bamFiles[@]}"

  samtools merge -f "${OUPUTDIR}/${GROUP}.merged.bam" "${bamFiles[@]}"
  samtools sort -@ 4 -o "${OUPUTDIR}/${GROUP}.sorted.merged.bam" "${OUPUTDIR}/${GROUP}.merged.bam"
  echo "Merged sample is sorted"
  samtools view -H "${OUPUTDIR}/${GROUP}.sorted.merged.bam" > "${OUPUTDIR}/${GROUP}.sorted.merged_header.sam"
  echo "Merged sample with samples that belong to ${GROUP} has been created."
  rm "${OUPUTDIR}/${GROUP}.merged.bam"
  
fi

## option PEAKS: peak calling is performed using MACS3 in paired-end at group level  
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
  
cat <<EOF

|| Running STEP 8.2 of PIPPA: PEAKS ||

Peak calling on ${GROUP} samples will be performed.
Output directory is ${PEAK_DIR_GROUPS}

EOF

  mkdir -p ${PEAK_DIR_GROUPS}

  sh "${SUB_SCRIPTS_DIR}/groupPeaks.sh" ${GROUP}

fi

## option COUNTS: get read counts in peaks called before
if [ $# = 1 ] || [[ $2 == 'COUNTS' ]]
then
  
  mkdir -p ${COUNTS_DIR_GROUPS}
  
cat <<EOF

|| Running STEP 8.3 of PIPPA: COUNTS ||

Counts in peaks will be output in a single file.
Differential accessible regions analysis will be performed to get topmost peaks.
Output directory is ${COUNTS_DIR_GROUPS}

EOF
 Rscript "${RSCRIPTS_DIR}/countsInPeaksGroup.r" ${CONFIGR} 

fi

## option CTCHECK: Results of previous stage are collated in an Rmarkdown  
if [ $# = 1 ] || [[ $2 == 'CTCHECK' ]]
then

cat <<EOF

|| Running STEP 8.4 of PIPPA: CTCHECK. Cell-type check results are collated in a report ||

Output directory will be: ${PEAKCOUNTS}

EOF
 
	Rscript -e "rmarkdown::render('${RSCRIPTS_DIR}/collateCellTypeCheck.Rmd', output_file='${PEAK_DIR}/QCOutput/stage3SummaryStats.html')" "${CONFIGR}" 
fi


conda deactivate

echo Job finished on:
date -u