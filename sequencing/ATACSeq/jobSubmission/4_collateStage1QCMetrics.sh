#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACQCSummaryStage1S4-%A.log
#SBATCH --error=ATACQCSummaryStage1S4-%A.err
#SBATCH --job-name=ATACQCSummaryStage1S4

## ===================================================================================================================##
##                          ATAC-seq pipeline STEP 4: Collate Stage 1 Results                                         ##
## ===================================================================================================================##
## EXECUTION: sbatch ./jobSubmission/4_collateStage1QCMetrics.sh <project-directory> <STEP>                           ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script collates results from the first part of the pipeline: pre-analysis and peak calling.      || 
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project dir> directory of project, location of config.txt file                                              ||
## $2 -> <STEP> Specify step to run: MULTIQC, COLLATE, SUMMARY. Can be combined. Default is to run all                ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, RAWDATADIR, ALIGNED_DIR, TRIM_DIR ||
##   SCRIPTS_DIR, PROJECT,FASTQCDIR                                                                                   ||
## - Version/directory of the following modules should be specified in config file: ACVERS, CONDA_ENV, CONDA, RVERS   ||
## - Softwares: Pandoc                                                                                                ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                            ||
## - Subscripts: progressReport.sh,  countMTcollateFS.sh                                                              ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R scripts: collateDataQualityStats.Rmd                                                                           ||
## - No array job number is required.                                                                                 ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time when script runs
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv ATACQCSummaryStage1S4-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR 
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir: $LOG_DIR

EOF


## ================ ##
##    VARIABLES     ##
## ================ ##

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "MULTIQC" ]] && [[ ! $2 =~ "COLLATE" ]] && [[ ! $2 =~ "SUMMARY" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use MULTIQC, COLLATE, SUMMARY or some combination of this as a single string (i.e. MULTIQC,COLLATE)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## option MULTIQC: collate QC output statistics in a single report   
if [ $# = 1 ] || [[ $2 =~ 'MULTIQC' ]]
then 
  
  module purge
  module load MultiQC

  mkdir -p ${FASTQCDIR}/multiqc
  mkdir -p ${ALIGNED_DIR}/multiqc
  mkdir -p ${TRIM_DIR}/multiqc
  
cat <<EOF

|| Running STEP 4.1 of ATAC-seq pipeline: MULTIQC ||

Results of fastqc on raw reads, trimming and alignment will be collated in a single report for each process.

Output dir for FASTQC multiqc is ${FASTQCDIR}/multiqc
Output dir for Alignment multiqc is ${ALIGNED_DIR}/multiqc
Output dir for Trimming multiqc is ${TRIM_DIR}/fastqc/multiqc

EOF
  
	cd ${FASTQCDIR}/
	multiqc . -f -o ${FASTQCDIR}/multiqc
	rm -f *.html
	
  cd ${ALIGNED_DIR}
	multiqc . -f -o ${ALIGNED_DIR}/multiqc
 
  cd ${TRIM_DIR}/qc
	multiqc . -f -o ${TRIM_DIR}/multiqc 
  rm -f ${TRIM_DIR}/qc/*.html
  
fi

## option COLLATE: collate other QC metrics
if [ $# = 1 ] || [[ $2 =~ 'COLLATE' ]]
then

cat <<EOF

|| Running step 4.2 of ATAC-seq pipeline: COLLATE ||

It will be checked that all samples had gone through all processes and had produced the right outputs.
Output directory is ${META_DIR} and ${ALIGNED_DIR}/ENCODEMetrics/

EOF
 
	sh ${SUB_SCRIPTS_DIR}/progressReport.sh 
	sh ${SUB_SCRIPTS_DIR}/countMTcollateFS.sh
fi

## option SUMMARY: collate output from previous steps into a single r markdown report
if [ $# = 1 ] || [[ $2 =~ 'SUMMARY' ]]
then
  
  module purge
	module load ${RVERS}
	module load Pandoc
 
cat <<EOF

|| Running step 4.3 of ATAC-seq pipeline: SUMMARY ||

Summary of results to this point are collated in a Rmarkdown report
Output directory is $PEAK_DIR/QCOutput

EOF
  
	Rscript -e "rmarkdown::render(paste0(commandArgs(trailingOnly=TRUE)[1],'/collateDataQualityStats.Rmd'), output_file=paste0(commandArgs(trailingOnly=TRUE)[2], '/QCOutput/stage1SummaryStats.html'))" "${RSCRIPTS_DIR}" "$PEAK_DIR" "${CONFIGR}"
fi

echo Job finished on:
date -u 
