#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/ATACQCSummary-%A.o
#SBATCH --error=ATACSeq/logFiles/%u/ATACQCSummary-%A.e
#SBATCH --job-name=ATACQCSummary

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder


## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR

#-----------------------------------------------------------------------#

if [ $# = 1 ] || [[ $2 =~ 'MULTIQC' ]]
then 
	module load MultiQC
	## use multiqc to collate QC output statistics

	mkdir -p ${FASTQCDIR}/multiqc
	cd ${FASTQCDIR}/
	multiqc . -f -o ${FASTQCDIR}/multiqc

	## remove redundant html files
	rm -f *.html
	rm -f ${TRIMDIR}/fastp_reports/*.html

	mkdir -p ${ALIGNEDDIR}/multiqc
	cd ${ALIGNEDDIR}/
	multiqc . -f -o ${ALIGNEDDIR}/multiqc
fi

if [ $# = 1 ] || [[ $2 =~ 'COLLATE' ]]
then
	## run other bespoke utilty scripts to collate other QC metrics
	cd ${SCRIPTDIR}

	./ATACSeq/preprocessing/progressReport.sh 
	./ATACSeq/preprocessing/countMTReads.sh 
	./ATACSeq/preprocessing/collateFlagStatOutput.sh 
fi


if [ $# = 1 ] || [[ $2 =~ 'SUMMARY' ]]
then
	## collate the earlier outputs into a r markdown report
	cd ${SCRIPTDIR}
	echo ${SCRIPTDIR}

	module load R/3.6.3-foss-2020a
	module load Pandoc
	Rscript -e "rmarkdown::render('ATACSeq/preprocessing/collateDataQualityStats.Rmd', output_file=paste0(commandArgs(trailingOnly=T)[1], '/QCOutput/stage1SummaryStats.html'))" "$PEAKDIR" "${SCRIPTDIR}" "$PROJECT" 
fi

shift #move command line arguments so that $1 is no longer project but step

if [[ $1 =~ 'FILTER' ]] #only run this if specified
then
	## collate the earlier outputs into a r markdown report
	cd ${SCRIPTDIR}

	shift

	module load R/3.6.3-foss-2020a
	Rscript ATACSeq/preprocessing/filterDataQualityStats.r ${PROJECT} $@ #all remaining cmd line arguments
fi

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
