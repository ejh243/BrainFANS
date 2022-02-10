#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/ATACAlignment-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment-%A_%a.e

## print start date and time
echo Job started on:
date -u

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

    step=$3
    all=1 #set to 1 to ensure if step flag is blank all steps are run
else
    step=$2
fi

## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}


## check step method matches required options and exit if not
if [[ ! $step =~ "FASTQC" ]] && [[ ! $step =~ "TRIM" ]] && [[ ! $step =~ "ALIGN" ]] && [[ ! $step =~ "ENCODE" ]] &&[[ ! $step == '' ]];
then 
    { echo "Unknown step specified. Please use FASTQC, TRIM, ALIGN, ENCODE or some combination of this as a single string (i.e. FASTQC,TRIM)" ; exit 1; }            
fi


## create array of all fastq files
cd ${RAWDATADIR}
FQFILES=($(find . -name '*[rR]1*q.gz')) ## this command searches for all fq files within

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""    

toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleID=$(basename ${toProcess%_[rR]*})
## later samples have an additional _S[num] in the file name need to remove
sampleID=${sampleID%_S[0-9]*}


## if number of flags is 1 (config.txt), then run all steps
if [ ${all} == 1 ] || [[ ${step} =~ 'FASTQC' ]]
then
    ## run sequencing QC and trimming on fastq files        
    module load FastQC 

    cd ${SCRIPTDIR}
    sh ./preScripts/fastqc.sh ${toProcess}  
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
    ## run alignment and post processing on sample
    module purge ## had conflict issues if this wasn't run first
    module load Bowtie2/2.3.4.2-foss-2018b
    module load SAMtools
    module load picard/2.6.0-Java-1.8.0_131

    cd ${SCRIPTDIR}
    sh ./ATACSeq/preprocessing/2_alignmentPE.sh ${toProcess}
fi

if [ ${all} == 1 ] || [[ ${step} =~ 'ENCODE' ]]
then
    module purge
    module load SAMtools
    module load picard/2.6.0-Java-1.8.0_131
    module load BEDTools/2.27.1-foss-2018b ##necessary to specify earlier BEDTools version to avoid conflict
    export PATH="$PATH:/gpfs/mrc0/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"

    ## load conda env for samstats
    module load Anaconda3
    source activate encodeqc

    cd ${SCRIPTDIR}
    sh ./ATACSeq/preprocessing/3_calcENCODEQCMetricsPE.sh ${sampleID}_sorted_chr1.bam
fi


## move log files into a folder
mkdir -p ATACSeq/logFiles/${USER}
mv ATACSeq/logFiles/ATACAlignment-${SLURM_ARRAY_JOB_ID}* ATACSeq/logFiles/${USER}
