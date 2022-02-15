#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/QCImputation.err # error file
#PBS -o LogFiles/QCImputation.log # output file


## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


cd $PBS_O_WORKDIR

####### 

## NOTE: Do not store confidenial information in this file use the config file

######

source ./SNPdata/config.txt
IMPUTATION=${DATADIR}/SNPdata/Merged/ImputationOutput

cd ${IMPUTATION}
module load R
Rscript ../../../scripts/SNPdata/summarizeImputation.r All/ ${KGG}/1000GP_Phase3_combined.legend ALL
Rscript ../../../scripts/SNPdata/summarizeImputation.r EUR/ ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab AF

module purge
module load VCFtools

sh combinedImputationOutput.sh

## print end date and time
echo Job finished:
date -u
