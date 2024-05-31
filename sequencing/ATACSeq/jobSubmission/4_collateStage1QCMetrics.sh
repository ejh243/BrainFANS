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
## EXECUTION: sbatch --array=0 ./sequencing/ATACSeq/jobSubmission/4_collateStage1QCMetrics.sh <project dir> <option>  ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array= can be any number                                                                                         ||
## $1 -> <project dir> project directory.                                                                             ||
## $2 -> <option> Specify step to run: MULTIQC, COLLATE, SUMMARY, COMPARE. Can be combined. Default is to run all     ||
##                                                                                                                    ||
## DESCRIPTION: This script collates results from the first part of the pipeline: pre-analysis and peak calling.      || 
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

## load config file provided on command line related to the specified project
source "${1}/config.txt"
echo "Loading config file for project: " $1
echo "Project directory is: " $DATADIR

LOG_DIR=ATACSeq/logFiles/${USER}/${SLURM_ARRAY_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACQCSummaryStage1S4-${SLURM_ARRAY_JOB_ID}* $LOG_DIR


## check step method matches required options and exit if not
if [[ ! $2 =~ "MULTIQC" ]] && [[ ! $2 =~ "COLLATE" ]] && [[ ! $2 =~ "SUMMARY" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use MULTIQC, COLLATE, SUMMARY or some combination of this as a single string (i.e. COMPARE,BATCH)" ; exit 1; }            
fi

## ============ ##
##    STEPS     ##
## ============ ##

## option MULTIQC: collate QC output statistics in a single report   
if [ $# = 1 ] || [[ $2 =~ 'MULTIQC' ]]
then 
  
  echo "Step 4.1 MULTIQC started. QC metrics will be collated"
	module load MultiQC

	mkdir -p ${FASTQCDIR}/multiqc
	cd ${FASTQCDIR}/
	multiqc . -f -o ${FASTQCDIR}/multiqc
 
	## remove redundant html files
	rm -f *.html
	rm -f ${TRIMDIR}/fastp_reports/*.html

	mkdir -p ${ALIGNEDDIR}/multiqc
  cd ${ALIGNEDDIR}
	multiqc . -f -o ${ALIGNEDDIR}/multiqc
 
  mkdir -p ${TRIMDIR}/multiqc
  cd ${TRIMDIR}/qc
	multiqc . -f -o ${TRIMDIR}/multiqc 
fi

## option COLLATE: collate other QC metrics
if [ $# = 1 ] || [[ $2 =~ 'COLLATE' ]]
then

  echo "Step 4.2 COLLATE started. QC metrics will be collated"
	
  cd ${SCRIPTDIR}
	sh ./ATACSeq/preprocessing/progressReport.sh 
	sh ./ATACSeq/preprocessing/countMTReads.sh 
	sh ./ATACSeq/preprocessing/collateFlagStatOutput.sh 
fi

## option SUMMARY: collate output from previous steps into a single r markdown report
if [ $# = 1 ] || [[ $2 =~ 'SUMMARY' ]]
then
  
  module purge
	module load ${RVERS}
	module load Pandoc
 
  echo "Step 4.3 SUMMARY started. Summary of results to this point are collated in a Rmarkdown report"
  cd ${SCRIPTDIR}
	Rscript -e "rmarkdown::render('ATACSeq/preprocessing/collateDataQualityStats.Rmd', output_file=paste0(commandArgs(trailingOnly=TRUE)[1], '/QCOutput/stage1SummaryStats.html'))" "$PEAKDIR" "${CONFIGR}"
fi

## option COMPARE: collate peak calling results into a single r markdown report
if [ $# = 1 ] || [[ $2 =~ 'COMPARE' ]]
then
  echo "Step 4.4 COMPARE started. Results from peak calling are compared and collated in two Rmarkdown reports: comparePeaksStats_1 and comparePeaksStats_2"
	cd ${SCRIPTDIR}
  
  module purge
	module load ${RVERS}
	module load Pandoc
	Rscript -e "rmarkdown::render('ATACSeq/preprocessing/collateComparePeaks_1.Rmd', output_file=paste0(commandArgs(trailingOnly=TRUE)[1], '/QCOutput/comparePeaksStats_1.html'))" "$PEAKDIR" "${CONFIGR}"
  Rscript -e "rmarkdown::render('ATACSeq/preprocessing/collateComparePeaks_2.Rmd', output_file=paste0(commandArgs(trailingOnly=TRUE)[1], '/QCOutput/comparePeaksStats_2.html'))" "$PEAKDIR" "${CONFIGR}"
fi

## option COMPARE: collate peak calling results into a single r markdown report
if [ $# = 1 ] || [[ $2 =~ 'BATCH' ]]
then
  echo "Step 4.4 COMPARE started. Results from peak calling are compared and collated in a Rmarkdown report"
	cd ${SCRIPTDIR}
  
  module purge
	module load ${RVERS}
	module load Pandoc
	Rscript -e "rmarkdown::render('ATACSeq/preprocessing/compareBatches.Rmd', output_file=paste0(commandArgs(trailingOnly=TRUE)[1], '/QCOutput/compareBatches.html'))" "$PEAKDIR" "${CONFIGR}"
fi
