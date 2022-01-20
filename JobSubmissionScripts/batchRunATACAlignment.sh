#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/ATAC/ATACAlignment-%A_%a.o
#SBATCH --error=LogFiles/ATAC/ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment-%A_%a.e
#SBATCH --array=0-410%50 

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

echo "Changing Folder to Data directory "
echo ${DATADIRPE}

cd ${DATADIRPE}
## find all folders with fastq files in
RAWDATADIR=($(find . -name '01_raw_reads'))

## create array of all fastq files
FQFILES=()
for FOLDER in ${RAWDATADIR[@]}
do     
    if [ ${FOLDER} != FastqParts ]
    then 
        echo "Processing folder: "
        echo ${FOLDER}
        FQFILES+=($(find ${FOLDER} -name '*[rR]1*q.gz')) ## this command searches for all fq files within
    fi
done

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""    


toProcess=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleID=$(basename ${toProcess%_[rR]*})
## later samples have an additional _S[num] in the file name need to remove
sampleID=${sampleID%_S[0-9]*}

## run sequencing QC and trimming on fastq files        
module load FastQC 
module load fastp

cd ${SCRIPTDIR}
sh ./ATACSeq/qcRawData.sh ${toProcess}  

## run alignment and post processing on sample
module purge ## had conflict issues if this wasn't run first
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131

cd ${SCRIPTDIR}
sh ./ATACSeq/alignmentPE.sh ${toProcess}

module purge
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
module load BEDTools/2.27.1-foss-2018b ##necessary to specify earlier BEDTools version to avoid conflict
export PATH="$PATH:/gpfs/mrc0/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"
export PATH="$PATH:/gpfs/mrc0/projects/Research_Project-MRC190311/software/SAMstats/bin/"


cd ${SCRIPTDIR}
sh ./ATACSeq/calcENCODEQCMetricsPE.sh ${sampleID}_sorted_chr1.bam