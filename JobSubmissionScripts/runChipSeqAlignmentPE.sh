#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=24:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/ChipAlignmentPE.err # error file
#PBS -o LogFiles/ChipAlignmentPE.log # output file

## needs to be executed from the scripts folder

## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

cd $PBS_O_WORKDIR
#cd ${SCRIPTDIR}
source ./ChipSeq/config.txt

## merge QC output
module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
cd ${DATADIRPE}
multiqc . -f ## can add flag to ignore certain folders if needed

## run alignment
module purge ## had conflict issues if this wasn't run first
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
#module load Java
sh ${SCRIPTDIR}/ChipSeq/alignmentPE.sh

echo Starting peak calling at:
date -u
module purge
module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
sh ${SCRIPTDIR}/ChipSeq/peakCallingPE.sh

## generate QC metrics
module purge
module load R
Rscript CalcChipQCMetrics.r 



## produce QC report


## print finish date and time
echo Job finished on:
date -u