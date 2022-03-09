#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=96:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=WGBS/logFiles/%u/WGBSAlignment-%A_%a.o
#SBATCH --error=WGBS/logFiles/%u/WGBSAlignment-%A_%a.e
#SBATCH --job-name=WGBSAlignment-%A_%a.e

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

source ./WGBS/config/config.txt 
echo "Project directory is: " $DATADIR


## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

echo "Changing Folder to Data directory "
echo ${RAWDATADIR}

cd ${RAWDATADIR}


## if raw data directory is empty, download files from ENA with specified ftp list
if [ -z "$(ls -A ${RAWDATADIR})" ]; 
then
   echo "Downloading files"
   ACCLIST=${METADIR}/file*
   python ${SCRIPTDIR}/preprocessing/ftp_url.py $ACCLIST #generate ena download list
   cd ${SCRIPTDIR}
   sh ./WGBS/preprocessing/_getENA.sh  
fi


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

      cd ${SCRIPTDIR}
      sh ./preScripts/fastqc.sh ${sampleID} ${toProcess[0]} ${toProcess[1]} 
	fi

	if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
	then
		module purge
		module load Trim_Galore

		cd ${SCRIPTDIR}
		echo "Changing to script directory: " ${SCRIPTDIR} ##
		sh ./preScripts/trimGalore.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}  

		echo "Finished Trim Galore on: " ##
		echo ${sampleID} ##
	fi

	if [ $# == 1 ] || [[ $2 =~ 'ALIGN' ]]
	then
		module purge
		module load Bismark

		cd ${SCRIPTDIR}
		sh ./WGBS/preprocessing/1_alignment.sh ${sampleID}
	fi

	echo 'EXITCODE: ' $?
	
	mkdir -p WGBS/logFiles/${USER}
	mv WGBS/logFiles/WGBSalign-${SLURM_ARRAY_JOB_ID}* WGBS/logFiles/${USER}
	
else
	echo "File list not found"
fi



## move log files into a folder
cd ${SCRIPTDIR}/WGBS/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv ATACAlignment-${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}

