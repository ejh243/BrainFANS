#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=DNAhydroxy/logFiles/CEGX5hmCQC-%A_%a.e # error file
#SBATCH --output=DNAhydroxy/logFiles/CEGX5hmCQC-%A_%a.o # output file
#SBATCH --job-name=CEGX5hmCQC-%A_%a.e



## print start date and time
echo Job started on:
date -u

## load config parameters
source $1

BAMFILES=($(ls ${ALIGNEDDIR}/*PC*bam))

echo "Number of bamfiles files found:"" ""${#BAMFILES[@]}"""    

toProcess=${BAMFILES[${SLURM_ARRAY_TASK_ID}]}
sampleID=$(basename ${toProcess%PC_*})

echo Starting peak calling at:
date -u

module load MACS2
module load BEDTools
module load SAMtools
module load Miniconda2
source activate epic2

##  run peak calling on each sample and calulate FRIP
./DNAhydroxy/preprocessing/1_samplePeaks.sh ${sampleID}

module purge
module load SAMtools

## extract sex chromosome for sex check
./DNAhydroxy/preprocessing/2_extractSexChr.sh ${sampleID}


## run rnaseqc to get proportion of reads assigned to genic regions
module load Miniconda2
source activate rnaseqc

QCDIR=${ALIGNEDDIR}/QCOutput

ICBAM=$(ls ${ALIGNEDDIR}/${sampleID}*IC*bam)

./${RNASEQCDIR}/rnaseqc.v2.4.2.linux ${GENCODEGTF/%annotation.gtf/genes.gtf} ${toProcess} ${QCDIR} -s ${sampleID}_PC
./${RNASEQCDIR}/rnaseqc.v2.4.2.linux ${GENCODEGTF/%annotation.gtf/genes.gtf} ${ICBAM} ${QCDIR} -s ${sampleID}_IC


mkdir -p DNAhydroxy/logFiles/${SLURM_ARRAY_JOB_ID}
mv DNAhydroxy/logFiles/CEGX5hmCQC-${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* DNAhydroxy/logFiles/${SLURM_ARRAY_JOB_ID}


## print finish date and time
echo Job finished on:
date -u