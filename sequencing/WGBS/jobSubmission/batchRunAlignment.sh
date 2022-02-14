#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=96:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=WGBS/logFiles/WGBSalign-%A_%a.o
#SBATCH --error=WGBS/logFiles/WGBSalign-%A_%a.e
#SBATCH --job-name=WGBSalign-%A_%a.e


## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "1. Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "2. Loading config file: "
echo $1
source ./$1 
all=$#

## if working in the development branch, load specified config.dev file
if [[ $2 =~ 'config.dev' ]]
then
    echo "Loading development config file:  "
    echo $2
    source ./$2

    step=$3
    all=1 #set to 1 to ensure if step flag is blank all steps are run
else
    step=$2
fi

## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $step =~ "FASTQC" ]] && [[ ! $step =~ "TRIM" ]] && [[ ! $step =~ "ALIGN" ]] && [[ ! $step == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

echo "3. Changing Folder to Data directory "
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
	
	if [ ${all} == 1 ] || [[ ${step} =~ 'FASTQC' ]]
	then
		## run sequencing QC on fastq files		
		module load FastQC/0.11.5-Java-1.7.0_80
		module load MultiQC
	
		cd ${SCRIPTDIR}
		echo "8. Changing to script directory: " ${SCRIPTDIR} ##
		sh ./WGBS/preprocessing/1_fastqc.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}
		echo "9. Finished fastqc on: " ##
		echo ${sampleID} ##
	fi

	if [ ${all} == 1 ] || [[ ${step} =~ 'TRIM' ]]
	then
		module purge
		module load Trim_Galore

		cd ${SCRIPTDIR}
		echo "8. Changing to script directory: " ${SCRIPTDIR} ##
		sh ./preScripts/trimGalore.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}  

		echo "9. Finished Trim Galore on: " ##
		echo ${sampleID} ##
	fi

	if [ ${all} == 1 ] || [[ ${step} =~ 'ALIGN' ]]
	then
		module purge
		module load Bismark

		cd ${SCRIPTDIR}
		sh ./WGBS/preprocessing/3_alignment.sh ${sampleID} ${toProcess[0]}   
	fi

	mkdir -p WGBS/logFiles/${USER}
	mv WGBS/logFiles/WGBSalign-${SLURM_ARRAY_JOB_ID}* WGBS/logFiles/${USER}
	
else
	echo "File list not found"
fi

