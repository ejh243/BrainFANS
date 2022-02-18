#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=ChIPSeq/logFiles/%u/ChIPAlignment-%A_%a.e # error file
#SBATCH --output=ChIPSeq/logFiles/%u/ChIPAlignment-%A_%a.o # output file
#SBATCH --job-name=ChIPAlignment

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


## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi


echo "Changing Folder to Data directory "
echo ${DATADIR}
cd ${DATADIR}

## find all R1 fastq files
FQFILES=()
FQFILES+=($(find . -name '*1.*q.gz')) ## this command searches for all fq files within ## rm [rR] add *_*1.*
echo ${FQFILES}

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

echo "SLURM_ARRAY_TASK_ID is: " "${SLURM_ARRAY_TASK_ID}"

toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${toProcess%.*.fastq.gz}) ##rm [rR]
echo "Current sample: " ${sampleName} 


if [ $# == 1 ] || [[ $2 =~ 'FASTQC' ]]
then
	## run sequencing QC on fastq files		
	module load FastQC/0.11.5-Java-1.7.0_80
	module load MultiQC
	module load fastp

	cd ${SCRIPTDIR}
	echo "Changing to script directory: " ${SCRIPTDIR} 
	sh ./preScripts/fastqc.sh ${toProcess}  

	echo "Finished fastqc on: " 
	echo ${sampleID} 
fi


if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
then
    module purge
    module load fastp
	
    cd ${SCRIPTDIR}
    echo "Changing to script directory: " ${SCRIPTDIR} ##
    sh ./preScripts/fastp.sh ${toProcess} 

    echo "Finished trim on: "
    echo ${sampleID}
fi


if [ $# == 1 ] || [[ $2 =~ 'ALIGN' ]]
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

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv ChIPAlignment-${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}/