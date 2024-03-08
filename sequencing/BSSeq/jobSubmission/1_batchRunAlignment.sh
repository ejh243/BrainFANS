#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=96:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=BSSeq/logFiles/%u/BSSeqAlignment-%A_%a.o
#SBATCH --error=BSSeq/logFiles/%u/BSSeqAlignment-%A_%a.e
#SBATCH --job-name=BSSeqAlignment-%A_%a.e

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for data type and project" $1
export PROJECT=$1

source ./BSSeq/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] && [[ ! $2 =~ "ENCODE" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

#-----------------------------------------------------------------------#

echo "Changing Folder to Data directory " ${RAWDATADIR}
cd ${RAWDATADIR}

## check if file containing list of sample IDs exists and if so:
if test -f ${METADIR}/samples.txt;
then 
	## create an array from the file
	mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 
	echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""

	sampleID=${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}
	## find the file name in RAWDATADIR
	toProcess=($(find ${RAWDATADIR} -maxdepth 1 -name ${SAMPLEIDS[${SLURM_ARRAY_TASK_ID}]}'*'))
	
	## sort the toProcess array so that R1 and R2 are consecutive 
	IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
	toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
	unset IFS 

	echo "R1 file found is: " $( basename ${toProcess[0]} )
	echo ${toProcess[0]}

	echo "Current sample: " ${sampleID} ##
	
	if [ $# == 1 ] || [[ $2 =~ 'FASTQC' ]]
	then
		## run sequencing QC on fastq files		     
      module load FastQC 

      mkdir -p ${FASTQCDIR}

      cd ${SCRIPTDIR}
      sh ./preScripts/fastqc.sh ${sampleID} ${toProcess[0]} ${toProcess[1]} 
	fi

	if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
	then
		module purge
		module load Trim_Galore

		mkdir -p ${TRIMDIR}

		cd ${SCRIPTDIR}
		echo "Changing to script directory: " ${SCRIPTDIR} ##
		sh ./preScripts/trimGalore.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}  

	fi

	if [ $# == 1 ] || [[ $2 =~ 'ALIGN' ]]
	then
		module purge
		module load Bismark

		mkdir -p $METHYLDIR

		cd ${SCRIPTDIR}
		sh ./BSSeq/preprocessing/alignment.sh ${sampleID}
	fi

	if [ $# == 1 ] || [[ $2 =~ 'ENCODE' ]]
	then
		module purge
		module load SAMtools
		module load BEDTools

		mkdir -p ${ALIGNEDDIR}/ENCODEMetrics

		cd ${SCRIPTDIR}
		sh ./BSSeq/preprocessing/calcQCMetrics.sh ${sampleID}
	fi

	echo 'EXITCODE: ' $?
	
	## move log files into a folder
	cd ${SCRIPTDIR}/BSSeq/logFiles/${USER}
	mkdir -p ${SLURM_ARRAY_JOB_ID}
	mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}
	
else
	echo "Sample list not found"
fi
