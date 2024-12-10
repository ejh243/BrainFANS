#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACcollateStage3QCS8-%A.log
#SBATCH --error=ATACcollateStage3QCS8-%A.err
#SBATCH --job-name=ATACcollateStage3QCS8

## ===================================================================================================================##
##                           ATAC-seq pipeline STEP 8: Cell type check                                                ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/8_collateStage3QCMetrics.sh <project directory> <STEP>                          ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the third stage of the ATAC-seq QC pipeline. The counts in peaks obtained in the ||
## previous step of the pipeline are used to confirm the cell-type identity of the samples. First, counts are         ||
## normalised and analysed (NORM). Then a report collating results and highlighting samples that have failed is       ||
## produced (CTCHECK).                                                                                                  ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: NORM, CHECK. Can be combined. Default is to run all                              ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR, PEAKCOUNTS ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA                                                       ||
## - A list of cell types of samples specified in config file: CELLTYPES                                              ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R subscripts: normCounts.r, collateCellTypeCheck.Rmd                                                             ||
## - A value for K parameter needs to be added to the R config file in order to run RUV-III-NB normalisation          ||
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
mv ATACcollateStage3QCS8-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

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

## check step method matches required options and exit if not
if [[ ! $2 =~ "NORM" ]] && [[ ! $2 =~ "CTCHECK" ]]  && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use NORM, CTCHECK, or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

## Input in command line what set of peaks to perform analysis on: promoter peaks (prom) or all peaks (all)
if [[ ! $3 == "PROM" ]] && [[ ! $3 == "ALL" ]];
then
  echo "No peak set specified. Please choose either PROM or ALL to select counts in promoter peaks or all peaks, respectively" ; exit 1; }     
fi

SETPEAKS=$3

## ========= ##
##   STEPS   ##
## ========= ##

## option NORM: Variance Partition Analysis of raw counts, normalisation of raw counts and again Variance partition analysis in normalised counts
if [ $# = 1 ] || [[ $2 == 'NORM' ]]
then

  mkdir -p ${PEAKCOUNTS}/normCounts/
  mkdir -p ${PEAKCOUNTS}/variance/
  
cat <<EOF

|| Running STEP 8.1 of ATAC-seq pipeline: NORM ||

Counts in peaks obtained in STEP 7.3 will be analysed and normalised. Variance Partition Analysis before and after normalisation is performed.

Output directory is ${PEAKCOUNTS}/variance
Output directory is ${PEAKCOUNTS}/normCounts

EOF

  Rscript "${RSCRIPTS_DIR}/normCounts.r" ${CONFIGR} ${SETPEAKS}
  
fi

## option GENCHECK: Sex check and Genotype check results are collated in Rmarkdown   
if [ $# = 1 ] || [[ $2 =~ 'CTCHECK' ]]
then

cat <<EOF

|| Running STEP 8.2 of ATAC-seq pipeline: CTCHECK. Results variance partition analysis, normalisation of counts and cell-type check are collated in a report ||

Output directory will be: ${PEAKCOUNTS}

EOF
 
	Rscript -e "rmarkdown::render(paste0(commandArgs(trailingOnly=TRUE)[1], '/collateCellTypeCheck.Rmd'), output_file=paste0(commandArgs(trailingOnly=TRUE)[2], '/stage3SummaryStats_all.html'))" "${RSCRIPTS_DIR}" "${PEAKCOUNTS}" "${CONFIGR}" "${SETPEAKS}"
  
fi
conda deactivate

echo Job finished on:
date -u