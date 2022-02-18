#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/ATAC/ATACAlignment-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/ATAC/ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment-%A_%a.e

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

## create array of all fastq files
cd ${RAWDATADIR}
FQFILES=($(find . -name '*[rR]1*q.gz')) ## this command searches for all fq files within

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""    

toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleID=$(basename ${toProcess%_[rR]*})
## later samples have an additional _S[num] in the file name need to remove
sampleID=${sampleID%_S[0-9]*}

if [ $# = 1 ] || [[ $2 =~ 'FASTQC' ]]
then
    ## run sequencing QC and trimming on fastq files        
    module load FastQC 

    cd ${SCRIPTDIR}
    sh ./preScripts/fastqc.sh ${toProcess}  
fi

if [ $# = 1 ] || [[ $2 =~ 'TRIM' ]]
then
    module purge
    module load fastp
	
    cd ${SCRIPTDIR}
    sh ./preScripts/fastp.sh ${toProcess} 
fi

if [ $# = 1 ] || [[ $2 =~ 'ALIGN' ]]
then
    ## run alignment and post processing on sample
    module purge ## had conflict issues if this wasn't run first
    module load Bowtie2/2.3.4.2-foss-2018b
    module load SAMtools
    module load picard/2.6.0-Java-1.8.0_131

    cd ${SCRIPTDIR}
    sh ./ATACSeq/preprocessing/2_alignment.sh ${toProcess}
fi

if [ $# = 1 ] || [[ $2 =~ 'ENCODE' ]]
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
    sh ./ATACSeq/preprocessing/3_calcENCODEQCMetrics.sh ${sampleID}_sorted_chr1.bam
fi

## move log files into a folder
mkdir -p logFiles/${SLURM_ARRAY_JOB_ID}
mv logFiles/ATACAlignment-${SLURM_ARRAY_JOB_ID}* logFiles/${SLURM_ARRAY_JOB_ID}
