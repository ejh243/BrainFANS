#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPAPeakSetCurationS9-%A.log
#SBATCH --error=PIPPAPeakSetCurationS9-%A.err
#SBATCH --job-name=PIPPAPeakSetCurationS9

## ===================================================================================================================##
##                                       PIPPA STEP 9: Curation of peak sets                                          ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/9_batchPeakSetCuration.sh <project directory> <STEP> <GROUP> <THRES>            ||
##                                                                                                                    ||
## DESCRIPTION: This script performs analysis and curation of peak sets that can be used in downstream analyses.      ||
## Samples that passed all QC are merged into single bam file (MERGE) and converted to BigWig format (BW). Peak       ||
## calling at cell type level is performed at default threshold (PEAKS). Segmented analysis is performed to guide     ||
## peak calling threshold (CUTOFF) for peak calling to be repeated (PEAKSTR).                                         ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: PEAKS. Can be combined. Default is to run all                                    ||
## $3 -> <GROUP> Specify cell group to run previous step(s).                                                          ||
## $4 -> <THRES> Specify peak calling threshold for curated peak set.                                                 ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt, config.r files in <project directory>.                                                               ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR_C,          ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA, CURATED_PEAKS, PEAK_DIR_CUTOFF                       ||
## - A list of cell types of samples specified in config files: CELLTYPES                                             ||
## - A conda environment setup with several modules: samtools, MACS3, R                                               ||
## - For other tools or references required, please refer to each subscript run in this script.                       ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: peaksForCuration.sh, peaksStrThres.sh                                                                ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: cutoffAnalysis.r, samplesForDownAnalysis.r                                                         ||
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
mv PIPPAPeakSetCurationS9-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

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

## check whether step to be run is specified
if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "MERGE" ]] && [[ ! $2 == "BW" ]] && [[ ! $2 == "PEAKS" ]] && [[ ! $2 == 'PEAKSTR' ]] && [[ ! $2 == 'CUTOFF' ]];
then 
    { echo "Unknown step specified. Please use PEAKS, BW or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## check cell type group has been specified for specific steps
if [[ ! " ${CELLTYPES[*]} " == *" $GROUP "* && ($2 == "MERGE" || $2 == "BW" || $2 == "PEAKS" || $2 == "CUTOFF" || $2 == "PEAKSTR" )]]; 
then
  { echo "Cell group specified not found. Please check whether samples that belong to that cell fraction exist."; exit 1; }            
fi

## check peak calling threshold has been specified for specific steps
if [[ $2 == "PEAKSTR" && -z "${THRES:-}" ]];
then 
    { echo "Invalid threshold for $2 STEP. Should be a peak calling threshold to be used to generate a curated peak set." ; exit 1; }
fi

if [[ ! ${GROUP} == '' ]]
then
  echo "Chosen cell fraction is ${GROUP}"
  ## Samples for group peak calling are gathered for input cell fraction
  #Rscript ${RSCRIPTS_DIR}/samplesForDownAnalysis.r ${CONFIGR} ${GROUP}

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

## option MERGE: merge bam files from samples in a single file
if [ $# == 1 ] || [[ $2 == 'MERGE' ]]
then

  OUPUTDIR="${ALIGNED_DIR}/merged/passAllQC/"
  mkdir -p "${OUPUTDIR}"
  
cat <<EOF

|| Running STEP 9.1 of ATAC-seq pipeline: MERGE ||

A merged bam file will be created with all samples for the chosen cell-type.
Output directory is ${OUPUTDIR}
Aligned files taken from ${ALIGNED_DIR}

EOF
  bamFiles=()
  for sampleName in ${SAMPLES[@]}
  do
    sampleBAM=${ALIGNED_DIR}/${sampleName}.filt.nodup.bam
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

## option BW: convert merged bam file to bigWig format
if [ $# == 1 ] || [[ $2 == 'BW' ]]
then
  
  INPUTDIR="${ALIGNED_DIR}/merged/passAllQC/"
  OUPUTDIR="${ALIGNED_DIR}/bigWig/passAllQC/"
  
  mkdir -p "${OUPUTDIR}"
  
cat <<EOF

|| Running STEP 9.2 of ATAC-seq pipeline: BW ||

If exists, a merged bam file will be converted to bigwig format to be used for genome tracks.
Output directory is ${OUPUTDIR}

EOF

  if test -f "${INPUTDIR}/${GROUP}.sorted.merged.bam";
  then
    samtools index "${INPUTDIR}/${GROUP}.sorted.merged.bam"
    
    bamCoverage \
    --bam "${INPUTDIR}/${GROUP}.sorted.merged.bam" \
    -o "${OUPUTDIR}/${GROUP}.bigWig" \
    --binSize 50 \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --minMappingQuality 30
  else
    echo "Merged bam file for ${GROUP} not found. Please refer to STEP 9.1"
  fi
fi

## option PEAKS: peak calling is performed using MACS3 in paired-end at group level in narrow peak mode 
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
  
  mkdir -p "${PEAK_DIR_C}/${GROUP}/"
  
cat <<EOF

|| Running STEP 9.3 of ATAC-seq pipeline: PEAKS ||

Peak calling on ${GROUP} samples will be performed.
Output directory is "${PEAK_DIR_C}/${GROUP}/"

EOF

  sh "${SUB_SCRIPTS_DIR}/peaksForCuration.sh" ${GROUP}
  
fi

## option CUTOFF: segmented analysis is performed on peak set to decide stringent threshold for curated peak set
if [ $# = 1 ] || [[ $2 == 'CUTOFF' ]]
then
  
  mkdir -p "${PEAK_DIR_CUTOFF}"
  
cat <<EOF

|| Running STEP 9.4 of ATAC-seq pipeline: CUTOFF ||

Segmented analysis is performed on peaks to decide peak calling threshold
Output directory is ${PEAK_DIR_CUTOFF}

EOF

  Rscript "${RSCRIPTS_DIR}/cutoffAnalysis.r" ${CONFIGR} ${GROUP}
  
fi

## option PEAKSTR: peak calling is repeated with chosen threshold for curated peak set
if [ $# = 1 ] || [[ $2 == 'PEAKSTR' ]]
then
  
  mkdir -p "${CURATED_PEAKS}"
  
cat <<EOF

|| Running STEP 9.5 of PIPPA: PEAKSTR ||

Peak calling on ${GROUP} samples will be performed.
Q value Cutoff is ${THRES}
Output directory is ${CURATED_PEAKS}

EOF

  sh "${SUB_SCRIPTS_DIR}/peaksStrThres.sh" ${GROUP} ${THRES}

fi

conda deactivate

echo Job finished on:
date -u