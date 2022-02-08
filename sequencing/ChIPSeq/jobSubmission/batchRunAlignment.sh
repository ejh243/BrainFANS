#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=ChIPSeq/logFiles/ChipAlignmentPE-%A_%a.e # error file
#SBATCH --output=ChIPSeq/logFiles/ChipAlignmentPE-%A_%a.o # output file
#SBATCH --job-name=ChipAlignmentPE

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1
all=$#

## if working in the development branch, load specified config.dev file
if [[ $2 =~ 'config.dev' ]]
then
    echo "Loading development config file:  "
    echo $2
    source ./$2

    step=$3 #as the last flag will be the steps to run flag
    all=1 #set to 1 to ensure if step flag is blank all steps are run
else
    step=$2
fi

## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $step =~ "FASTQC" ]] && [[ ! $step =~ "TRIM" ]] && [[ ! $step =~ "ALIGN" ]] &&[[ ! $step == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi


echo "Changing Folder to Data directory "
echo ${DATADIRPE}
cd ${DATADIRPE}

## find all R1 fastq files
FQFILES=()
FQFILES+=($(find . -name '*1.*q.gz')) ## this command searches for all fq files within ## rm [rR] add *_*1.*
echo ${FQFILES}

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

echo "SLURM_ARRAY_TASK_ID is: " "${SLURM_ARRAY_TASK_ID}"

toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${toProcess%.*.fastq.gz}) ##rm [rR]
echo "Current sample: " ${sampleName} ##


if [ ${all} == 1 ] || [[ ${step} =~ 'FASTQC' ]]
then
	## run sequencing QC on fastq files		
	module load FastQC/0.11.5-Java-1.7.0_80
	module load MultiQC
	module load fastp

	cd ${SCRIPTDIR}
	echo "8. Changing to script directory: " ${SCRIPTDIR} ##
	sh ./preScripts/fastqc.sh ${toProcess}  

	echo "9. Finished fastqc on: " ##
	echo ${sampleID} ##
fi


if [ ${all} == 1 ] || [[ ${step} =~ 'TRIM' ]]
then
    module purge
    module load fastp
	
    cd ${SCRIPTDIR}
    sh ./preScripts/fastp.sh ${toProcess} 
fi


if [ ${all} == 1 ] || [[ ${step} =~ 'ALIGN' ]]
then
	module purge ## had conflict issues if this wasn't run first
	module load Bowtie2
	module load SAMtools
	module load picard/2.6.0-Java-1.8.0_131
	module load BEDTools
	module load Java
	cd ${SCRIPTDIR}
	sh ./ChIPSeq/preprocessing/2_alignment.sh ${toProcess} ## using ./ rather than sh executes script in current session and can make use of variables alredy declared.
fi

