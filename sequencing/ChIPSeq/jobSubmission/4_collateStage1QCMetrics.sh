#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ChIPSeq/logFiles/%u/ChIPQCSummary-%A.o
#SBATCH --error=ChIPSeq/logFiles/%u/ChIPQCSummary-%A.e
#SBATCH --job-name=ChIPQCSummary

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ChIPSeq/config/config.txt 
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

	sh ./ChIPSeq/preprocessing/progressReport.sh 
	sh ./ChIPSeq/preprocessing/collateFlagStatOutput.sh 
fi

if [ $# = 1 ] || [[ $2 =~ 'SUMMARY' ]]
then
	## collate the earlier outputs into a r markdown report
	cd ${SCRIPTDIR}
	echo ${SCRIPTDIR}

	module load R/3.6.3-foss-2020a
	module load Pandoc
	Rscript -e "rmarkdown::render('ChIPSeq/preprocessing/collateS1SumStats.Rmd', output_file=paste0(commandArgs(trailingOnly=T)[1], '/QCOutput/stage1SummaryStats.html'))" "$PEAKDIR" "${SCRIPTDIR}" "$PROJECT" 
fi

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_JOB_ID}
mv *${SLURM_JOB_ID}.* ${SLURM_JOB_ID}
