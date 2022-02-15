#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/QCImputation.err # error file
#PBS -o LogFiles/QCImputation.log # output file



## print start date and time
echo Job started on:
date -u


####### 

## NOTE: Do not store confidenial information in this file use the config file

######

source ./$1


cd ${IMPUTEDIR}/Output

module load R/3.5.1-foss-2018b-Python-2.7.15

Rscript SNPArray/preprocessing/5_summarizeImputation.r All/ ${KGG}/1000GP_Phase3_combined.legend ALL
Rscript SNPArray/preprocessing/5_summarizeImputation.r EUR/ ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab AF

module purge
module load VCFtools

## combine imputation output separately for ALL and EUR versions
sh SNPArray/preprocessing/6_combineImputationOutput.sh ${IMPUTEDIR}/Output/All
sh SNPArray/preprocessing/6_combineImputationOutput.sh ${IMPUTEDIR}/Output/EUR

## print end date and time
echo Job finished:
date -u
