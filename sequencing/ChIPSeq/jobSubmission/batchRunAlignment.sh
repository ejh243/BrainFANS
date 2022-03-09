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
JOBNAME="ChIPAlignment"
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ChIPSeq/config/config.txt 
echo "Project directory is: " $DATADIR


## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] && [[ ! $2 =~ "ENCODE" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi


echo "Changing Folder to Data directory "
echo ${DATADIR}
cd ${DATADIR}

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

   ## if number of flags is 1 (config.txt), then run all steps
    if [ $# == 1 ] || [[ $2 =~ 'FASTQC' ]]
    then
        ## run sequencing QC and trimming on fastq files        
        module load FastQC 

        cd ${SCRIPTDIR}
        sh ./preScripts/fastqc.sh ${sampleID} ${toProcess[0]} ${toProcess[1]}  
    fi

    if [ $# == 1 ] || [[ $2 =~ 'TRIM' ]]
    then
        module purge
        module load fastp
    	
        cd ${SCRIPTDIR}
        sh ./preScripts/fastp.sh ${sampleID} ${toProcess[0]} ${toProcess[1]} 
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
		sh ./ChIPSeq/preprocessing/2_alignment.sh ${sampleID} ## using ./ rather than sh executes script in current session and can make use of variables alredy declared.
	fi

	if [ $# = 1 ] || [[ $2 =~ 'ENCODE' ]]
	then
		module purge
        ## load conda env for samstats
        module load Anaconda3
        source activate encodeqc
        module load SAMtools
        module load BEDTools/2.27.1-foss-2018b ##necessary to specify earlier BEDTools version
        module load Java

        cd ${SCRIPTDIR}
        sh ./ChIPSeq/preprocessing/3_calcENCODEQCMetrics.sh ${sampleID}
    fi

	echo 'EXITCODE: ' $?

	## move log files into a folder
	cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
	mkdir -p ${SLURM_ARRAY_JOB_ID}
	mv ${JOBNAME}-${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}/
else
	echo 'File list not found'
fi