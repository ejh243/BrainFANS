#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACGroupPeakCallingS7-%A-%a.o
#SBATCH --error=ATACGroupPeakCallingS7-%A-%a.e
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
## $2 -> <STEP> Specify step to run: PEAKS, FILT, COUNTS, DIFFCOUNTS and CTCHECK. Can be combined.                    ||
##  Default is to run all                                                                                             ||
## $3 -> <GROUP> Specify cell group to run FILT or PEAKS steps. Not needed for other steps.                           ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR_GROUPS,     ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA, PEAKCOUNTS_DIR_GROUPS                                ||
## - A list of cell types of samples specified in config file: CELLTYPES                                              ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, MACS3                                                  ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: groupPeaks.sh                                                                                        ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: countsInPeaksGroup.r, samplesForGroupAnalysis.r ,diffPeakGroupsCounts.r ,collateCellTypeCheck.Rmd  ||
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
if [[ ! $2 =~ "PEAKS" ]] && [[ ! $2 =~ "COUNTS" ]] && [[ ! $2 == 'FILT' ]] && [[ ! $2 == 'CTCHECK' ]] && [[ ! $2 == 'DIFFCOUNTS' ]] && [[ ! $2 == ' ' ]];
then 
    { echo "Unknown step specified. Please use FRAGSIZE, PEAKS, COUNTS or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

if [[ ! " ${CELLTYPES[*]} " == *" $GROUP "* && ($2 =~ "FRAGSIZE" || $2 =~ "PEAKS" || $2 =~ "FILT")]]; 
then
  { echo "Cell group specified not found. Please check whether samples that belong to that cell fraction exist."; exit 1; }            
fi


if [[ ! ${GROUP} == '' && ( $2 =~ "PEAKS" || $2 =~ "FILT")]];
then
  echo "Chosen cell fraction is ${GROUP}"
  ## Samples for group peak calling are gathered for input cell fraction
  Rscript ${RSCRIPTS_DIR}/samplesForGroupAnalysis.r ${CONFIGR} ${GROUP}

  ## Retrieve samples from the same cell fraction
  mapfile -t SAMPLES < ${META_DIR}/samplesForGroupAnalysisOrdered_${GROUP}.txt

else
  echo "No cell group specified."
fi


## ========= ##
##   STEPS   ##
## ========= ##

## option PEAKS: peak calling is performed using MACS3 in paired-end at group level  
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
  
cat <<EOF

|| Running STEP 7.1 of ATAC-seq pipeline: PEAKS ||

Peak calling on ${GROUP} samples will be performed.
Output directory is ${PEAK_DIR_GROUPS}

EOF

  mkdir -p ${PEAK_DIR_GROUPS}

  sh "${SUB_SCRIPTS_DIR}/groupPeaks.sh" ${GROUP}

fi

## option FILT: peaks called at cell-fraction level are filtered based on peaks called at sample level
if [ $# = 1 ] || [[ $2 == 'FILT' ]]
then
  
cat <<EOF

|| Running STEP 7.2 of ATAC-seq pipeline: FILT ||

Filtering of ${GROUP} peaks will be performed.
Output directory is ${PEAK_DIR_GROUPS}

EOF

  Rscript "${RSCRIPTS_DIR}/filtGroupPeaks.r" ${CONFIGR} ${GROUP}

fi

## option COUNTS: get read counts in peaks called above
if [ $# = 1 ] || [[ $2 == 'COUNTS' ]]
then
  
  mkdir -p ${PEAKCOUNTS_DIR_GROUPS}
  
cat <<EOF

|| Running STEP 7.3 of ATAC-seq pipeline: COUNTS ||

Counts in peaks will be output in a single file.
Output directory is ${PEAKCOUNTS_DIR_GROUPS}

EOF
  Rscript "${RSCRIPTS_DIR}/countsInPeaksGroup.r" ${CONFIGR} 
fi

## option DIFFCOUNTS: get top most differential peaks between cell-types
if [ $# = 1 ] || [[ $2 == 'DIFFCOUNTS' ]]
then
  
  mkdir -p ${PEAKCOUNTS_DIR_GROUPS}
  
cat <<EOF

|| Running STEP 7.4 of ATAC-seq pipeline: DIFFCOUNTS ||

Differential analysis between cell-types will be performed.
Previous STEP 7.3 should be run before.

Output directory is ${PEAKCOUNTS_DIR_GROUPS}

EOF

  Rscript "${RSCRIPTS_DIR}/diffPeakGroupsCounts.r" ${CONFIGR} 
fi

## option CTCHECK: Results of previous stage are collated in an Rmarkdown  
if [ $# = 1 ] || [[ $2 == 'CTCHECK' ]]
then

cat <<EOF

|| Running STEP 7.5 of ATAC-seq pipeline: CTCHECK. Results variance partition analysis, normalisation of counts and cell-type check are collated in a report ||

Output directory will be: ${PEAKCOUNTS}

EOF
 
	Rscript -e "rmarkdown::render('${RSCRIPTS_DIR}/collateCellTypeCheck.Rmd', output_file='${PEAKCOUNTS}/stage3SummaryStats.html')" "${CONFIGR}" 
fi

conda deactivate

echo Job finished on:
date -u
