#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=BSSeq/logFiles/%u/BSSeqQCSummary-%A.o
#SBATCH --error=BSSeq/logFiles/%u/BSSeqQCSummary-%A.e
#SBATCH --job-name=BSSeqQCSummary


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

source ./BSSeq/config/config.txt 
echo "Project directory is: " $DATADIR


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
	echo 'tbc'
	cd ${SCRIPTDIR}

	sh ./BSSeq/preprocessing/collateENCODEQCOutput.sh
	sh ./BSSeq/preprocessing/progressReport.sh 
fi

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/BSSeq/logFiles/${USER}
mkdir -p ${SLURM_JOB_ID}
mv BSSeqQCSummary-${SLURM_JOB_ID}* ${SLURM_JOB_ID}/
