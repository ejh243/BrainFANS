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
echo ${DATADIR}

cd ${DATADIR}
## find all folders with fastq files in

## if raw data directory is empty, download files from ENA with specified ftp list
if [ -z "$(ls -A ${DATADIR})" ]; 
then
   echo "Downloading files"
   ACCLIST=${METADIR}/file*
   python ${SCRIPTDIR}preprocessing/ftp_url.py $ACCLIST #generate ena download list
   cd ${SCRIPTDIR}
   sh ./WGBS/preprocessing/_getENA.sh  
fi

## create array of all fastq files
cd ${DATADIR}
FQFILES=()
FQFILES+=($(find . -name '*[_rR]1*q.gz')) ## this command searches for all fq files within 
echo ${FQFILES}


echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

echo "SLURM_ARRAY_TASK_ID is: " "${SLURM_ARRAY_TASK_ID}"

toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleID=$(basename ${toProcess%_*}) ##rm [rR]
echo "7. Current sample: " ${sampleID} ##

if [ ${all} == 1 ] || [[ ${step} =~ 'FASTQC' ]]
then
	## run sequencing QC and trimming on fastq files		
	module load FastQC/0.11.5-Java-1.7.0_80
	module load MultiQC
	##module load fastp

	cd ${SCRIPTDIR}
	echo "8. Changing to script directory: " ${SCRIPTDIR} ##
	sh ./WGBS/preprocessing/1_fastqc.sh ${toProcess}  
	echo "9. Finished fastqc on: " ##
	echo ${sampleID} ##
fi

if [ ${all} == 1 ] || [[ ${step} =~ 'TRIM' ]]
then
	module purge
	module load Trim_Galore

	cd ${SCRIPTDIR}
	echo "8. Changing to script directory: " ${SCRIPTDIR} ##
	sh ./preScripts/trimGalore.sh ${toProcess}  

	echo "9. Finished Trim Galore on: " ##
	echo ${sampleID} ##
fi

if [ ${all} == 1 ] || [[ ${step} =~ 'ALIGN' ]]
then
	module purge
	module load Bismark

	cd ${SCRIPTDIR}
	sh ./WGBS/preprocessing/3_alignment.sh ${toProcess}
fi

mkdir -p WGBS/logFiles/WGBSalign_${SLURM_ARRAY_JOB_ID}
mv WGBS/logFiles/WGBSalign-${SLURM_ARRAY_JOB_ID}* WGBS/logFiles/WGBSalign_${SLURM_ARRAY_JOB_ID}


