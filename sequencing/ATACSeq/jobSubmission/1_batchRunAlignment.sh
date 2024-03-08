#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/ATACAlignment-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
    
## needs to be executed from the scripts folder

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
source ./ATACSeq/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 =~ "FASTQC" ]] && [[ ! $2 =~ "TRIM" ]] && [[ ! $2 =~ "ALIGN" ]] && [[ ! $2 =~ "ENCODE" ]] &&[[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi

#-----------------------------------------------------------------------#

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


    ## if number of flags is 1 ($PROJECT), then run all steps
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
		## run alignment and post processing on sample
		module purge ## had conflict issues if this wasn't run first
		module load Bowtie2/2.3.4.2-foss-2018b
		module load SAMtools
		module load picard/2.6.0-Java-1.8.0_131
		export PATH="$PATH:/lustre/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"
		
		cd ${SCRIPTDIR}
		sh ./ATACSeq/preprocessing/alignment.sh ${sampleID}
	fi

	if [ $# == 1 ] || [[ $2 =~ 'ENCODE' ]]
	then
		module purge
		module load SAMtools
		module load BEDTools/2.27.1-foss-2018b ##necessary to specify earlier BEDTools version
        
        ## load conda env for samstats
        module load Anaconda3
        source activate encodeqc

        cd ${SCRIPTDIR}
        sh ./ATACSeq/preprocessing/calcENCODEQCMetrics.sh ${sampleID}
    fi

    ## move log files into a folder
    cd ATACSeq/logFiles/${USER}/
    mkdir -p ${SLURM_ARRAY_JOB_ID}
    mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}


else
    echo "File list not found"
fi

echo 'EXIT CODE: ' $?

